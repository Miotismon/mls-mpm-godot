using Godot;
using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;

public partial class MLSMPM3DFluidMultithreadGPU : Node3D
{
    [StructLayout(LayoutKind.Sequential)]
    struct Particle
    {
        public Vector3 pos; // position
        public float padding_pos; // padding for uniform
        public Vector3 vel; // velocity
        public float mass;
        public Vector3 C_x; // affine matrix
        public float padding_C_x;
        public Vector3 C_y; 
        public float padding_C_y;
        public Vector3 C_z; 
        public float padding_C_z;
    }
    int particle_struct_size = Marshal.SizeOf<Particle>();


    [StructLayout(LayoutKind.Sequential)]
    struct Cell // with fixed point numbers
    {
        public int vel_x;
        public int vel_y;
        public int vel_z;
        public int mass;
    }
    int cell_struct_size = Marshal.SizeOf<Cell>();

    //[StructLayout(LayoutKind.Sequential)]
    //struct ParticlePos
    //{
    //    public Vector3 pos;
    //    public float padding;
    //}
    //int particle_pos_struct_size = Marshal.SizeOf<ParticlePos>();

    Vector3I grid_size = new Vector3I(64, 64, 64);  // TODO: make grid_size actually a vec3 in the shaders, rn it's just a int that takes grid_size.x for all dimensions
    int num_cells;

    const int max_particle_count = 300000;
    int num_particles;

    Particle[] ps;
    Cell[] grid;


    // simulation parameters
    [ExportGroup("Simulation Parameters")]

    // dt = the time step of our simulation
    private float dt = 0.2f;
    [Export(PropertyHint.Range, "0.0f,0.4f,")]
    float Dt
    {
        get;
        set
        {
            field = Math.Clamp(value, 0.0f, 0.4f);
            UpdatePushConstants();
        }
    } = 0.2f;
    [Export]
    int sim_iterations = 2;
    [Export]
    float gravity = -0.3f;

    // fluid parameters
    //[ExportSubgroup("Fluid Parameters")]
    [Export]
    float rest_density = 4.0f;
    [Export]
    float dynamic_viscosity = 0.1f;
    
    // equation of state
    [Export]
    float eos_stiffness = 1.0f;
    [Export]
    float eos_power = 7.0f;


    [ExportGroup("Setup")]
    // interaction
    [Export]
    PhysicsBody3D sphere_body;
    bool mouse_down = false;
    Vector2 mouse_pos_local;
    Vector3 sphere_pos;
    float sphere_radius = 15.0f;
    Vector3 mouse_on_plane_pos;

    // fixed point numbers
    static int fixed_point_mult = (int)1e7;

    // rendering
    MultiMeshInstance3D multi_mesh_instance;
    //GpuParticles3D gpu_particles;
    //MeshInstance3D screen_space_quad;

    // compute shaders
    [Export(PropertyHint.File, "*.glsl")]
    String clear_grid_shader_file_path;
    [Export(PropertyHint.File, "*.glsl")]
    String p2g_1_shader_file_path;
    [Export(PropertyHint.File, "*.glsl")]
    String p2g_2_shader_file_path;
    [Export(PropertyHint.File, "*.glsl")]
    String update_grid_shader_file_path;
    [Export(PropertyHint.File, "*.glsl")]
    String g2p_shader_file_path;

    const int workgroup_size = 128;
    uint num_grid_workgroups;
    uint num_particle_workgroups;

    uint particle_pos_tex_width;
    uint particle_pos_tex_height;

    RenderingDevice rd;

    Rid grid_buffer;
    Rid particle_buffer;
    Rid particle_pos_tex;

    Rid clear_grid_shader_rid;
    //byte[] clear_grid_push_constants;
    Rid clear_grid_uniform_set;
    Rid clear_grid_pipeline;

    Rid p2g_1_shader_rid;
    byte[] p2g_1_push_constants;
    Rid p2g_1_uniform_set;
    Rid p2g_1_pipeline;

    Rid p2g_2_shader_rid;
    byte[] p2g_2_push_constants;
    Rid p2g_2_uniform_set;
    Rid p2g_2_pipeline;

    Rid update_grid_shader_rid;
    byte[] update_grid_push_constants;
    Rid update_grid_uniform_set;
    Rid update_grid_pipeline;

    Rid g2p_shader_rid;
    byte[] g2p_push_constants;
    Rid g2p_uniform_set;
    Rid g2p_pipeline;


    

    public override void _Ready()
    {
        // onready variables
        multi_mesh_instance = GetNode<MultiMeshInstance3D>("MultiMeshInstance3D");
        //gpu_particles = GetNode<GpuParticles3D>("GPUParticles3D");
        //screen_space_quad = GetNode<MeshInstance3D>("FullScreenQuad");

        if (sphere_body != null)
        {
            sphere_pos = sphere_body.GlobalPosition;
        } 
        else
        {
            GD.PrintErr("sphere_body not set");
        }

        num_cells = grid_size.X * grid_size.Y * grid_size.Z;
        

        
        InitialiseSim();

        // prepare for particles
        if (multi_mesh_instance != null)
        {
            multi_mesh_instance.Multimesh.InstanceCount = num_particles;
            for (int i = 0; i < num_particles; i++)
            {
                multi_mesh_instance.Multimesh.SetInstanceTransform(i, Transform3D.Identity);
            }
        }
        
        //if (gpu_particles != null)
        //{
        //    gpu_particles.Amount = num_particles;
        //}
        

        particle_pos_tex_width = (uint)Mathf.Sqrt(num_particles) + 1;
        particle_pos_tex_height = particle_pos_tex_width;


        InitGPU();

        // set globals for rendering etc.
        var global_node = GetTree().Root.GetNode<Node>("Global");
        global_node.Set("particle_count", num_particles);
        global_node.Set("particle_pos_texture", particle_pos_tex);
        global_node.Set("particle_pos_texture_width", particle_pos_tex_width);
        global_node.Set("current_simulator", this);


        //var outputBytes = rd.BufferGetData(particle_buffer);
        //var i = 0;
        //Vector3 particle_pos = new Vector3(
        //    BitConverter.ToSingle(outputBytes, i * particle_struct_size + 0),
        //    BitConverter.ToSingle(outputBytes, i * particle_struct_size + 4),
        //    BitConverter.ToSingle(outputBytes, i * particle_struct_size + 8)
        //);

        //GD.Print("particle pos after: ", particle_pos);



        //// Read back the data from the buffers
        //var outputBytes = rd.BufferGetData(particle_buffer);
        //var output = new float[ps.Length * Marshal.SizeOf<Particle>()];
        //Buffer.BlockCopy(outputBytes, 0, output, 0, outputBytes.Length);

        //GD.Print($"Particle 1 pos_xyz: {output[0]},{output[1]},{output[2]}");
        //GD.Print($"at gridpos 1,0,10, old value: {grid[10].vel_y}; grid_res: {grid_res}; new value: {output[((1 * grid_res * grid_res) + (0 * grid_res) + 10) * 4 + 1]} ");

    }



    public override void _Process(double delta)
    {
        
        HandleMouseInteraction();

        //GD.Print($"process delta: {delta}");
        //GD.Print($"process delta: {delta}  dt: {Dt},  sim_iterations: {sim_iterations}");
        for (int i = 0; i < sim_iterations; ++i)
        {
            SetComputeLists();

            // only use if using local rendering device
            //rd.Submit();
            //rd.Sync();

        }

    }

    public void TurnOnVisualisation()
    {
        multi_mesh_instance.Visible = true;
    }

    public void TurnOffVisualisation()
    {
        multi_mesh_instance.Visible = false;
    }


    // These resources are expensive to make, so create them once and cache for subsequent runs
    void InitGPU()
    {

        // get global rendering device, but make sure not to delete it when this node gets deleted
        rd = RenderingServer.GetRenderingDevice();
        if (rd == null)
        {
            OS.Alert($"Couldn't get global RenderingDevice on GPU: {RenderingServer.GetVideoAdapterName()}\n \nNote: RenderingDevice is only available in the Forward + and Mobile rendering methods, not Compatibility.");
            return;
        }



        // load glsl shaders
        if (clear_grid_shader_file_path == "")
        {
            GD.PrintErr("No Shader File export given");
            return;
        }
        clear_grid_shader_rid = LoadShader(ref rd, clear_grid_shader_file_path);
        p2g_1_shader_rid = LoadShader(ref rd, p2g_1_shader_file_path);
        p2g_2_shader_rid = LoadShader(ref rd, p2g_2_shader_file_path);
        update_grid_shader_rid = LoadShader(ref rd, update_grid_shader_file_path);
        g2p_shader_rid = LoadShader(ref rd, g2p_shader_file_path);



        // copy data to byte array
        int particle_array_size = particle_struct_size * ps.Length;

        byte[] particle_bytes = new byte[particle_array_size];
        IntPtr ptr = Marshal.AllocHGlobal(particle_array_size);

        try
        {
            // Copy struct array to unmanaged memory
            for (int i = 0; i < ps.Length; i++)
            {
                IntPtr dest = IntPtr.Add(ptr, i * particle_struct_size);
                Marshal.StructureToPtr(ps[i], dest, false);
            }

            // Copy unmanaged memory to byte array
            Marshal.Copy(ptr, particle_bytes, 0, particle_array_size);
        }
        finally
        {
            Marshal.FreeHGlobal(ptr);
        }



        // create storage buffers
        int grid_array_size = cell_struct_size * grid.Length;
        grid_buffer = rd.StorageBufferCreate((uint)grid_array_size);

        //int particle_array_size = particle_struct_size * num_particles;
        particle_buffer = rd.StorageBufferCreate((uint)particle_bytes.Length, particle_bytes);

        //int particle_pos_size = particle_pos_struct_size * num_particles;
        //particle_pos_buffer = rd.StorageBufferCreate((uint)particle_pos_size);

        //Vector3 debug_particle_pos = new Vector3(
        //    BitConverter.ToSingle(particle_bytes, 0),
        //    BitConverter.ToSingle(particle_bytes, 4),
        //    BitConverter.ToSingle(particle_bytes, 8)
        //);
        //Vector3 debug_particle_vel = new Vector3(
        //    BitConverter.ToSingle(particle_bytes, 16),
        //    BitConverter.ToSingle(particle_bytes, 20),
        //    BitConverter.ToSingle(particle_bytes, 24)
        //);
        //float debug_particle_mass = BitConverter.ToSingle(particle_bytes, 28);
        //GD.Print($"Particle data in Buffer at start: pos = {debug_particle_pos}, vel = {debug_particle_vel}, mass = {debug_particle_mass}");


        // create texture format
        RDTextureFormat particle_pos_format = new RDTextureFormat
        {
            Format = RenderingDevice.DataFormat.R32G32B32A32Sfloat,
            TextureType = RenderingDevice.TextureType.Type2D,
            Width = particle_pos_tex_width,
            Height = particle_pos_tex_height,
            Depth = 1,
            Mipmaps = 1,
            ArrayLayers = 1,
            UsageBits = RenderingDevice.TextureUsageBits.StorageBit | RenderingDevice.TextureUsageBits.SamplingBit
        };

        // create texture
        particle_pos_tex = rd.TextureCreate(particle_pos_format, new RDTextureView());


        // create uniforms to assign the storage buffers and texture to the rendering device
        RDUniform grid_uniform = new RDUniform
        {
            UniformType = RenderingDevice.UniformType.StorageBuffer,
            Binding = 0
        };
        grid_uniform.AddId(grid_buffer);

        RDUniform particle_uniform = new RDUniform
        {
            UniformType = RenderingDevice.UniformType.StorageBuffer,
            Binding = 1
        };
        particle_uniform.AddId(particle_buffer);

        RDUniform particle_pos_uniform = new RDUniform
        {
            UniformType = RenderingDevice.UniformType.Image,
            Binding = 2
        };
        particle_pos_uniform.AddId(particle_pos_tex);


        UpdatePushConstants();

        // add uniforms to uniform sets for each shader
        clear_grid_uniform_set = rd.UniformSetCreate([grid_uniform], clear_grid_shader_rid, 0);
        p2g_1_uniform_set = rd.UniformSetCreate([grid_uniform, particle_uniform], p2g_1_shader_rid, 0);
        p2g_2_uniform_set = rd.UniformSetCreate([grid_uniform, particle_uniform], p2g_2_shader_rid, 0);
        update_grid_uniform_set = rd.UniformSetCreate([grid_uniform], update_grid_shader_rid, 0);
        g2p_uniform_set = rd.UniformSetCreate([grid_uniform, particle_uniform, particle_pos_uniform], g2p_shader_rid, 0);

        // create compute pipeline for each shader
        clear_grid_pipeline = rd.ComputePipelineCreate(clear_grid_shader_rid);
        p2g_1_pipeline = rd.ComputePipelineCreate(p2g_1_shader_rid);
        p2g_2_pipeline = rd.ComputePipelineCreate(p2g_2_shader_rid);
        update_grid_pipeline = rd.ComputePipelineCreate(update_grid_shader_rid);
        g2p_pipeline = rd.ComputePipelineCreate(g2p_shader_rid);


        // configure workgroups
        num_grid_workgroups = ((uint)num_cells - 1 / workgroup_size) + 1;
        num_particle_workgroups = (uint)(num_particles - 1 / workgroup_size) + 1;


        // set references to particle position texture if needed
        if (multi_mesh_instance != null)
        {
            ShaderMaterial multi_mesh_instance_material = (ShaderMaterial)multi_mesh_instance.Multimesh.Mesh.SurfaceGetMaterial(0);
            // set uniforms in material
            multi_mesh_instance_material.SetShaderParameter("particle_pos_tex_width", particle_pos_tex_width);
            // replace reference to particle_pos_tex in material, only works if using the global rendering device
            Texture2Drd particle_pos_tex_rd0 = (Texture2Drd)multi_mesh_instance_material.GetShaderParameter("particle_pos_tex");
            particle_pos_tex_rd0.TextureRdRid = particle_pos_tex;
        }
        
        //if (gpu_particles != null)
        //{ 
        //    // get shader material from gpu particle
        //    ShaderMaterial gpu_particle_material = (ShaderMaterial)gpu_particles.ProcessMaterial;
        //    // set uniforms in material
        //    gpu_particle_material.SetShaderParameter("particle_pos_tex_width", particle_pos_tex_width);
        //    // replace reference to particle_pos_tex in material, only works if using the global rendering device
        //    Texture2Drd particle_pos_tex_rd1 = (Texture2Drd)gpu_particle_material.GetShaderParameter("particle_pos_tex");
        //    particle_pos_tex_rd1.TextureRdRid = particle_pos_tex;
        //}

        //if (screen_space_quad != null)
        //{
        //    ShaderMaterial screen_space_quad_material = (ShaderMaterial)screen_space_quad.GetSurfaceOverrideMaterial(0);
        //    screen_space_quad_material.SetShaderParameter("particle_pos_tex_width", particle_pos_tex_width);
        //    Texture2Drd particle_pos_tex_rd2 = (Texture2Drd)screen_space_quad_material.GetShaderParameter("particle_pos_tex");
        //    particle_pos_tex_rd2.TextureRdRid = particle_pos_tex;
        //}
        
        

    }

    Rid LoadShader(ref RenderingDevice rd, String path)
    {
        RDShaderFile shader_file_data = (RDShaderFile)GD.Load(path);
        RDShaderSpirV shader_spirv = shader_file_data.GetSpirV();
        return rd.ShaderCreateFromSpirV(shader_spirv);
    }

    private void UpdatePushConstants()
    {
        // setup push constants for each shader
        byte[] fixed_point_mult_bytes = BitConverter.GetBytes(fixed_point_mult);
        byte[] grid_res_bytes = BitConverter.GetBytes(grid_size.X);
        byte[] dt_bytes = BitConverter.GetBytes(Dt);
        byte[] gravity_bytes = BitConverter.GetBytes(gravity);

        byte[] rest_density_bytes = BitConverter.GetBytes(rest_density);
        byte[] dynamic_viscosity_bytes = BitConverter.GetBytes(dynamic_viscosity);
        byte[] eos_stiffness_bytes = BitConverter.GetBytes(eos_stiffness);
        byte[] eos_power_bytes = BitConverter.GetBytes(eos_power);

        byte[] sphere_pos_bytes = new byte[12];
        Buffer.BlockCopy(BitConverter.GetBytes(sphere_pos.X), 0, sphere_pos_bytes, 0, 4);
        Buffer.BlockCopy(BitConverter.GetBytes(sphere_pos.Y), 0, sphere_pos_bytes, 4, 4);
        Buffer.BlockCopy(BitConverter.GetBytes(sphere_pos.Z), 0, sphere_pos_bytes, 8, 4);
        byte[] tex_width_bytes = BitConverter.GetBytes(particle_pos_tex_width);

        if (p2g_1_push_constants == null)
        {
            p2g_1_push_constants = new byte[fixed_point_mult_bytes.Length + grid_res_bytes.Length + 8]; // add 8 bytes to round up to multiple of 16 for GLSL shader

        }
        Buffer.BlockCopy(fixed_point_mult_bytes, 0, p2g_1_push_constants, 0, fixed_point_mult_bytes.Length);
        Buffer.BlockCopy(grid_res_bytes, 0, p2g_1_push_constants, fixed_point_mult_bytes.Length, grid_res_bytes.Length);

        if (p2g_2_push_constants == null)
        {
            p2g_2_push_constants = new byte[fixed_point_mult_bytes.Length + grid_res_bytes.Length + dt_bytes.Length + rest_density_bytes.Length + dynamic_viscosity_bytes.Length + eos_stiffness_bytes.Length + eos_power_bytes.Length + 4];

        }
        Buffer.BlockCopy(fixed_point_mult_bytes, 0, p2g_2_push_constants, 0, fixed_point_mult_bytes.Length);
        Buffer.BlockCopy(grid_res_bytes, 0, p2g_2_push_constants, 4, grid_res_bytes.Length);
        Buffer.BlockCopy(dt_bytes, 0, p2g_2_push_constants, 8, dt_bytes.Length);
        Buffer.BlockCopy(rest_density_bytes, 0, p2g_2_push_constants, 12, rest_density_bytes.Length);
        Buffer.BlockCopy(dynamic_viscosity_bytes, 0, p2g_2_push_constants, 16, dynamic_viscosity_bytes.Length);
        Buffer.BlockCopy(eos_stiffness_bytes, 0, p2g_2_push_constants, 20, eos_stiffness_bytes.Length);
        Buffer.BlockCopy(eos_power_bytes, 0, p2g_2_push_constants, 24, eos_power_bytes.Length);

        if (update_grid_push_constants == null)
        {
            update_grid_push_constants = new byte[fixed_point_mult_bytes.Length + grid_res_bytes.Length + dt_bytes.Length + gravity_bytes.Length];
        }
        Buffer.BlockCopy(fixed_point_mult_bytes, 0, update_grid_push_constants, 0, fixed_point_mult_bytes.Length);
        Buffer.BlockCopy(grid_res_bytes, 0, update_grid_push_constants, 4, grid_res_bytes.Length);
        Buffer.BlockCopy(dt_bytes, 0, update_grid_push_constants, 8, dt_bytes.Length);
        Buffer.BlockCopy(gravity_bytes, 0, update_grid_push_constants, 12, gravity_bytes.Length);

        if(g2p_push_constants == null)
        {
            g2p_push_constants = new byte[fixed_point_mult_bytes.Length + grid_res_bytes.Length + dt_bytes.Length + 4 + sphere_pos_bytes.Length + tex_width_bytes.Length];

        }
        Buffer.BlockCopy(fixed_point_mult_bytes, 0, g2p_push_constants, 0, fixed_point_mult_bytes.Length);
        Buffer.BlockCopy(grid_res_bytes, 0, g2p_push_constants, 4, grid_res_bytes.Length);
        Buffer.BlockCopy(dt_bytes, 0, g2p_push_constants, 8, dt_bytes.Length);
        Buffer.BlockCopy(sphere_pos_bytes, 0, g2p_push_constants, 16, sphere_pos_bytes.Length);
        Buffer.BlockCopy(tex_width_bytes, 0, g2p_push_constants, 28, tex_width_bytes.Length);
    }

    void SetComputeLists()
    {
        //GD.Print("seting up compute lists");

        // setup compute lists
        var clear_grid_compute_list = rd.ComputeListBegin();
        rd.ComputeListBindComputePipeline(clear_grid_compute_list, clear_grid_pipeline);
        //rd.ComputeListSetPushConstant(clear_grid_compute_list, clear_grid_push_constants, (uint)clear_grid_push_constants.Length);
        rd.ComputeListBindUniformSet(clear_grid_compute_list, clear_grid_uniform_set, 0);
        rd.ComputeListDispatch(clear_grid_compute_list, num_grid_workgroups, 1, 1);
        rd.ComputeListEnd();

        var p2g_1_compute_list = rd.ComputeListBegin();
        rd.ComputeListBindComputePipeline(p2g_1_compute_list, p2g_1_pipeline);
        rd.ComputeListSetPushConstant(p2g_1_compute_list, p2g_1_push_constants, (uint)p2g_1_push_constants.Length);
        rd.ComputeListBindUniformSet(p2g_1_compute_list, p2g_1_uniform_set, 0);
        rd.ComputeListDispatch(p2g_1_compute_list, num_particle_workgroups, 1, 1);
        rd.ComputeListEnd();

        var p2g_2_compute_list = rd.ComputeListBegin();
        rd.ComputeListBindComputePipeline(p2g_2_compute_list, p2g_2_pipeline);
        rd.ComputeListSetPushConstant(p2g_2_compute_list, p2g_2_push_constants, (uint)p2g_2_push_constants.Length);
        rd.ComputeListBindUniformSet(p2g_2_compute_list, p2g_2_uniform_set, 0);
        rd.ComputeListDispatch(p2g_2_compute_list, num_particle_workgroups, 1, 1);
        rd.ComputeListEnd();

        var update_grid_compute_list = rd.ComputeListBegin();
        rd.ComputeListBindComputePipeline(update_grid_compute_list, update_grid_pipeline);
        rd.ComputeListSetPushConstant(update_grid_compute_list, update_grid_push_constants, (uint)update_grid_push_constants.Length);
        rd.ComputeListBindUniformSet(update_grid_compute_list, update_grid_uniform_set, 0);
        rd.ComputeListDispatch(update_grid_compute_list, num_grid_workgroups, 1, 1);
        rd.ComputeListEnd();

        var g2p_compute_list = rd.ComputeListBegin();
        rd.ComputeListBindComputePipeline(g2p_compute_list, g2p_pipeline);
        rd.ComputeListSetPushConstant(g2p_compute_list, g2p_push_constants, (uint)g2p_push_constants.Length);
        rd.ComputeListBindUniformSet(g2p_compute_list, g2p_uniform_set, 0);
        rd.ComputeListDispatch(g2p_compute_list, num_particle_workgroups, 1, 1);
        rd.ComputeListEnd();
    }

    void CleanupGpu()
    {
        //if (particle_pos_tex == new Rid())
        //{
        //    return;
        //}

        // all resources must be freed after use to avoid memory leaks.

        // pipelines
        rd.FreeRid(clear_grid_pipeline);
        clear_grid_pipeline = new Rid();

        rd.FreeRid(p2g_1_pipeline);
        p2g_1_pipeline = new Rid();

        rd.FreeRid(p2g_2_pipeline);
        p2g_2_pipeline = new Rid();

        rd.FreeRid(update_grid_pipeline);
        update_grid_pipeline = new Rid();

        rd.FreeRid(g2p_pipeline);
        g2p_pipeline = new Rid();

        // uniform sets
        rd.FreeRid(clear_grid_uniform_set);
        clear_grid_uniform_set = new Rid();

        rd.FreeRid(p2g_1_uniform_set);
        p2g_1_uniform_set = new Rid();

        rd.FreeRid(p2g_2_uniform_set);
        p2g_2_uniform_set = new Rid();

        rd.FreeRid(update_grid_uniform_set);
        update_grid_uniform_set = new Rid();

        rd.FreeRid(g2p_uniform_set);
        g2p_uniform_set = new Rid();

        // shader rids
        rd.FreeRid(clear_grid_shader_rid);
        clear_grid_shader_rid = new Rid();

        rd.FreeRid(p2g_1_shader_rid);
        p2g_1_shader_rid = new Rid();

        rd.FreeRid(p2g_2_shader_rid);
        p2g_2_shader_rid = new Rid();

        rd.FreeRid(update_grid_shader_rid);
        update_grid_shader_rid = new Rid();

        rd.FreeRid(g2p_shader_rid);
        g2p_shader_rid = new Rid();

        // buffers
        rd.FreeRid(grid_buffer);
        grid_buffer = new Rid();

        rd.FreeRid(particle_buffer);
        particle_buffer = new Rid();

        rd.FreeRid(particle_pos_tex);
        particle_pos_tex = new Rid();

        // rendering device, DON'T DO THIS IF YOU USE THE GLOBAL RENDERING DEVICE
        //rd.Free();
        //rd = null;
    }

    void HandleMouseInteraction()
    {
        mouse_down = false;
        if (Input.IsMouseButtonPressed(MouseButton.Left) && sphere_body != null)
        {
            mouse_down = true;
            //var mouse_pos_global = GetGlobalMousePosition();
            //mouse_pos_local = (new Vector2(mouse_pos_global.X, mouse_pos_global.Y)) - GlobalPosition;

            //var camera3D = GetNode<Camera3D>("Camera3D");
            //var from = camera3D.ProjectRayOrigin(eventMouseButton.Position);
            //var to = from + camera3D.ProjectRayNormal(eventMouseButton.Position) * RayLength;

            sphere_pos = mouse_on_plane_pos;
            sphere_body.GlobalPosition = sphere_pos;
            //GD.Print("sphere_pos: ", sphere_pos);

            byte[] sphere_pos_bytes = new byte[12];
            Buffer.BlockCopy(BitConverter.GetBytes(sphere_pos.X), 0, sphere_pos_bytes, 0, 4);
            Buffer.BlockCopy(BitConverter.GetBytes(sphere_pos.Y), 0, sphere_pos_bytes, 4, 4);
            Buffer.BlockCopy(BitConverter.GetBytes(sphere_pos.Z), 0, sphere_pos_bytes, 8, 4);

            Buffer.BlockCopy(sphere_pos_bytes, 0, g2p_push_constants, 16, sphere_pos_bytes.Length);
        }
    }

    public static int EncodeFixedPoint(float floating_point)
    {
        return (int)(floating_point * fixed_point_mult);
    }

    public static float DecodeFixedPoint(int fixed_point)
    {
        return (float)(fixed_point) / fixed_point_mult;
    }

    void InitialiseSim()
    {
        // initialising a bunch of points in a box
        List<Vector3> temp_positions = new List<Vector3>();
        const float spacing = 0.6f;
        const int box_x = 32, box_y = 32, box_z = 32;
        float sx = grid_size.X / 2.0f, sy = grid_size.Y / 2.0f, sz = grid_size.Z / 2.0f;
        for (float i = sx - box_x / 2; i < sx + box_x / 2; i += spacing)
        {
            for (float j = sy - box_y / 2; j < sy + box_y / 2; j += spacing)
            {
                for (float k = sz - box_z / 2; k < sz + box_z / 2; k += spacing)
                {
                    Vector3 pos = new Vector3(i, j, k);
                    temp_positions.Add(pos);
                }
            }
        }
        //GD.Print("temp_positions.Count: ", temp_positions.Count);

        // round number of particles to nearest power of 2
        //int po2_num = 1;
        //while (po2_num <= temp_positions.Count)
        //{
        //    po2_num <<= 1;
        //}
        //num_particles = po2_num >> 1;
        num_particles = temp_positions.Count;

        ps = new Particle[num_particles];

        // populate our array of particles, set their initial state
        for (int i = 0; i < num_particles; ++i)
        {
            Particle p = new Particle();
            p.pos = temp_positions[i];
            p.vel = Vector3.Zero;
            p.mass = 1.0f;
            p.C_x = Vector3.Zero;
            p.C_y = Vector3.Zero;
            p.C_z = Vector3.Zero;
            ps[i] = p;
        }

        grid = new Cell[num_cells];

        for (int i = 0; i < num_cells; ++i)
        {
            var cell = new Cell();
            grid[i] = cell;
        }

        GD.Print("Number of Particles: ", num_particles);
    }

    public override void _Notification(int what)
    {
        if (what == NotificationPredelete)
        {
            CleanupGpu();
        }
    }

    public void On_sphere_move_plane_input_event(Node camera, InputEvent new_event, Vector3 event_position, Vector3 normal, int shape_idx)
    {
        if (new_event is InputEventMouseMotion)
        {
            mouse_on_plane_pos = event_position;
            //GD.Print("event position: ", event_position);
        }
    }

}
