using Godot;
using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using System.Threading;
using System.Threading.Tasks;

public partial class MLSMPM3DFluidMultithreadNew : Node3D
{
    struct Particle
    {
        public Vector3 pos; // position
        public Vector3 vel; // velocity
        public Basis C; // affine matrix
        public float mass;
    }

    struct Cell // with fixed point numbers
    {
        public int vel_x;
        public int vel_y;
        public int vel_z;
        public int mass;
    }

    const int grid_res = 32;
    const int num_cells = grid_res * grid_res * grid_res;


    // simulation parameters


    // dt = the time step of our simulation
    const float dt = 0.2f;
    const float sim_iterations = (int)(1.0f / dt);

    const float gravity = -0.3f;

    // fluid parameters
    const float rest_density = 4.0f;
    const float dynamic_viscosity = 0.1f;
    // equation of state
    const float eos_stiffness = 10.0f;
    const float eos_power = 4;


    int num_particles;

    Particle[] ps;
    Cell[] grid;

    // fixed point numbers
    static int fixed_point_mult = (int)1e7;


    // interaction
    //const float mouse_radius = 10.0f;
    bool mouse_down = false;
    Vector2 mouse_pos_local;
    [Export]
    PhysicsBody3D sphere_body;
    Vector3 sphere_pos = Vector3.Zero;
    float sphere_radius = 15f;
    Vector3 mouse_on_plane_pos;

    // rendering stuff
    MultiMeshInstance3D multiMeshInstanceParticles;
    //MultiMeshInstance2D multiMeshInstanceCells;

    public override void _Ready()
    {
        // onready variables
        multiMeshInstanceParticles = GetNode<MultiMeshInstance3D>("MultiMeshInstanceParticles");
        //multiMeshInstanceCells = GetNode<MultiMeshInstance2D>("MultiMeshInstanceCells");
        if (sphere_body != null)
        {
            sphere_pos = sphere_body.GlobalPosition;
        }

        

        Initialise();

        // prepare multimesh for particles
        multiMeshInstanceParticles.Multimesh.InstanceCount = num_particles;

        // draw out all grid cells in multimesh
        //multiMeshInstanceCells.Multimesh.InstanceCount = grid.Length;
        //for (int x = 0; x < grid_res; x++)
        //{
        //    for (int y = 0; y < grid_res; y++)
        //    {
        //        QuadMesh quadMesh = new QuadMesh();
        //        quadMesh.Size = new Vector2(1.0f, -1.0f);
        //        multiMeshInstanceCells.Multimesh.Mesh = quadMesh;
        //        multiMeshInstanceCells.Multimesh.SetInstanceTransform2D(x * grid_res + y, new Transform2D(0.0f, new Vector2(x, y)));
        //    }
        //}

    }



    public override void _Process(double delta)
    {
        HandleMouseInteraction();

        //GD.Print("dt: " + dt + "; sim_iterations: " + sim_iterations);
        for (int i = 0; i < sim_iterations; ++i)
        {
            Simulate();
        }


        // update particle positions in multimesh
        for (int i = 0; i < num_particles; i++)
        {
            Particle particle = ps[i];

            //multiMeshInstanceParticles.Multimesh.SetInstanceTransform2D(i, new Transform2D(0.0f, particle.pos));
            multiMeshInstanceParticles.Multimesh.SetInstanceTransform(i, new Transform3D(Basis.Identity, particle.pos));
            //multiMeshInstanceParticles.Multimesh.SetInstanceTransform(i, new Transform3D());
        }

    }

    public override void _PhysicsProcess(double delta)
    {

    }

    void HandleMouseInteraction()
    {
        mouse_down = false;
        if (Input.IsMouseButtonPressed(MouseButton.Left))
        {
            mouse_down = true;
            //var mouse_pos_global = GetGlobalMousePosition();
            //mouse_pos_local = (new Vector2(mouse_pos_global.X, mouse_pos_global.Y)) - GlobalPosition;

            //var camera3D = GetNode<Camera3D>("Camera3D");
            //var from = camera3D.ProjectRayOrigin(eventMouseButton.Position);
            //var to = from + camera3D.ProjectRayNormal(eventMouseButton.Position) * RayLength;

            sphere_pos = mouse_on_plane_pos;
            sphere_body.GlobalPosition = sphere_pos;
            GD.Print("sphere_pos: ", sphere_pos);
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

    void Initialise()
    {
        // initialising a bunch of points in a box
        List<Vector3> temp_positions = new List<Vector3>();
        const float spacing = 1.0f;
        const int box_x = 16, box_y = 16, box_z = 16;
        const float sx = grid_res / 2.0f, sy = grid_res / 2.0f, sz = grid_res / 2.0f;
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

        
        num_particles = temp_positions.Count;

        ps = new Particle[num_particles];

        // populate our array of particles, set their initial state
        for (int i = 0; i < num_particles; ++i)
        {
            Particle p = new Particle();
            p.pos = temp_positions[i];
            p.vel = Vector3.Zero;
            p.C = new Basis();
            p.mass = 1.0f;
            ps[i] = p;
        }

        grid = new Cell[num_cells];

        for (int i = 0; i < num_cells; ++i)
        {
            var cell = new Cell();
            grid[i] = cell;
        }


    }

    void Simulate()
    {


        // reset grid scratchpad
        var start = Time.GetTicksUsec();
        ClearGrid();
        var end = Time.GetTicksUsec();
        var cleargrid_time = (end - start) / 1000.0;

        // P2G first round
        start = Time.GetTicksUsec();
        P2G_1();
        end = Time.GetTicksUsec();
        var p2g_1_time = (end - start) / 1000.0;

        // P2G, second round
        start = Time.GetTicksUsec();
        P2G_2();
        end = Time.GetTicksUsec();
        var p2g_2_time = (end - start) / 1000.0;

        // grid velocity update
        start = Time.GetTicksUsec();
        UpdateGrid();
        end = Time.GetTicksUsec();
        var update_grid_time = (end - start) / 1000.0;

        // G2P
        start = Time.GetTicksUsec();
        G2P();
        end = Time.GetTicksUsec();
        var g2p_time = (end - start) / 1000.0;

        //GD.Print(String.Format("ClearGrid: {0}ms, P2G_1: {1}ms, P2G_2: {2}ms, UpdateGrid: {3}, G2P: {4}ms, num_particles: {5}", cleargrid_time, p2g_1_time, p2g_2_time, update_grid_time, g2p_time, num_particles));
    }

    private void ClearGrid()
    {
        var options = new ParallelOptions
        {
            MaxDegreeOfParallelism = -1
        };

        Parallel.For(0, num_cells, options, i =>
        {
            ClearCell(i);
        });

        //Task[] tasks = new Task[num_cells];

        //for (int i = 0; i < num_cells; ++i)
        //{
        //    tasks.Append(Task.Run(i => ClearCell(i)));
        //}
    }

    private void ClearCell(int i)
    {
        var cell = grid[i];

        cell.mass = 0;
        cell.vel_x = 0;
        cell.vel_y = 0;
        cell.vel_z = 0;

        grid[i] = cell;
    }

    private void P2G_1()
    {
        var options = new ParallelOptions
        {
            MaxDegreeOfParallelism = -1
        };

        Parallel.For(0, num_particles, options, i =>
        {
            P2G_1_Single(i);
        });
    }

    private void P2G_1_Single(int i)
    {
        var p = ps[i];

        // quadratic interpolation weights
        Vector3I cell_idx = (Vector3I)p.pos;
        Vector3 cell_diff = (p.pos - cell_idx) - new Vector3(0.5f, 0.5f, 0.5f);
        Vector3[] weights = {
                0.5f * ((new Vector3(0.5f, 0.5f, 0.5f) - cell_diff) * (new Vector3(0.5f, 0.5f, 0.5f) - cell_diff)),
                new Vector3(0.75f, 0.75f, 0.75f) - (cell_diff * cell_diff),
                0.5f * ((new Vector3(0.5f, 0.5f, 0.5f) + cell_diff) * (new Vector3(0.5f, 0.5f, 0.5f) + cell_diff))
        };

        Basis C = p.C;

        for (uint gx = 0; gx < 3; ++gx)
        {
            for (uint gy = 0; gy < 3; ++gy)
            {
                for (uint gz = 0; gz < 3; ++gz)
                {
                    float weight = weights[gx].X * weights[gy].Y * weights[gz].Z;

                    Vector3I cell_x = new Vector3I((int)(cell_idx.X + gx - 1), (int)(cell_idx.Y + gy - 1), (int)(cell_idx.Z + gz - 1));
                    Vector3 cell_dist = (cell_x - p.pos) + new Vector3(0.5f, 0.5f, 0.5f);
                    Vector3 Q = C * cell_dist;
                    
                    float mass_contrib = weight * p.mass;
                    Vector3 vel_contrib = mass_contrib * (p.vel + Q);

                    // cell_index = x*width*depth + y*width + z
                    int cell_index = (int)cell_x.X * grid_res * grid_res + (int)cell_x.Y * grid_res + (int)cell_x.Z;

                    // mass and momentum update
                    //Cell cell = grid[cell_index];

                    //cell.mass += mass_contrib;
                    //cell.vel += vel_contrib;

                    //grid[cell_index] = cell;

                    //grid[cell_index].mass += EncodeFixedPoint(mass_contrib);
                    //grid[cell_index].vel_x += EncodeFixedPoint(vel_contrib.X);
                    //grid[cell_index].vel_y += EncodeFixedPoint(vel_contrib.Y);
                    //grid[cell_index].vel_z += EncodeFixedPoint(vel_contrib.Z);

                    Interlocked.Add(ref grid[cell_index].mass, EncodeFixedPoint(mass_contrib));
                    Interlocked.Add(ref grid[cell_index].vel_x, EncodeFixedPoint(vel_contrib.X));
                    Interlocked.Add(ref grid[cell_index].vel_y, EncodeFixedPoint(vel_contrib.Y));
                    Interlocked.Add(ref grid[cell_index].vel_z, EncodeFixedPoint(vel_contrib.Z));
                }
            }
        }
    }

    private void P2G_2()
    {
        var options = new ParallelOptions
        {
            MaxDegreeOfParallelism = -1
        };

        Parallel.For(0, num_particles, options, i =>
        {
            P2G_2_Single(i);
        });
    }

    private void P2G_2_Single(int i)
    {


        var p = ps[i];

        // quadratic interpolation weights
        Vector3I cell_idx = (Vector3I)p.pos;
        Vector3 cell_diff = (p.pos - cell_idx) - new Vector3(0.5f, 0.5f, 0.5f);
        Vector3[] weights = {
                0.5f * ((new Vector3(0.5f, 0.5f, 0.5f) - cell_diff) * (new Vector3(0.5f, 0.5f, 0.5f) - cell_diff)),
                new Vector3(0.75f, 0.75f, 0.75f) - (cell_diff * cell_diff),
                0.5f * ((new Vector3(0.5f, 0.5f, 0.5f) + cell_diff) * (new Vector3(0.5f, 0.5f, 0.5f) + cell_diff))
        };

        float density = 0.0f;
        uint gx, gy, gz;
        for (gx = 0; gx < 3; ++gx)
        {
            for (gy = 0; gy < 3; ++gy)
            {
                for (gz = 0; gz < 3; ++gz)
                {
                    float weight = weights[gx].X * weights[gy].Y * weights[gz].Z;
                    int cell_index = (int)(cell_idx.X + gx - 1) * grid_res * grid_res + (int)(cell_idx.Y + gy - 1) * grid_res + (int)(cell_idx.Z + gz - 1);
                    density += DecodeFixedPoint(grid[cell_index].mass) * weight;
                }
                    
            }
        }

        float volume = p.mass / density;

        float pressure = Mathf.Max(-0.1f, eos_stiffness * (Mathf.Pow(density / rest_density, eos_power) - 1));

        Basis stress = new Basis(
            new Vector3( -pressure, 0.0f,       0.0f    ),
            new Vector3( 0.0f,      -pressure,  0.0f    ),
            new Vector3( 0.0f,      0.0f,       -pressure)
        );

        // velocity gradient
        Basis dudv = p.C;
        Basis dudvT = dudv.Transposed();
        Basis strain = new Basis(dudv.X + dudvT.X, dudv.Y + dudvT.Y, dudv.Z + dudvT.Z);

        Basis viscosity_term = new Basis(strain.X * dynamic_viscosity, strain.Y * dynamic_viscosity, strain.Z * dynamic_viscosity);
        stress = new Basis(stress.X + viscosity_term.X, stress.Y + viscosity_term.Y, stress.Z + viscosity_term.Z);

        Basis eq_16_term_0 = new Basis(-volume * 4 * stress.X * dt, -volume * 4 * stress.Y * dt, -volume * 4 * stress.Z * dt);//-volume * 4 * stress * dt;

        for (gx = 0; gx < 3; ++gx)
        {
            for (gy = 0; gy < 3; ++gy)
            {
                for (gz = 0; gz < 3; ++gz)
                {
                    float weight = weights[gx].X * weights[gy].Y * weights[gz].Z;

                    Vector3I cell_x = new Vector3I((int)(cell_idx.X + gx - 1), (int)(cell_idx.Y + gy - 1), (int)(cell_idx.Z + gz - 1));
                    Vector3 cell_dist = (cell_x - p.pos) + new Vector3(0.5f, 0.5f, 0.5f);

                    int cell_index = (int)cell_x.X * grid_res * grid_res + (int)cell_x.Y * grid_res + (int)cell_x.Z;

                    //Cell cell = grid[cell_index];

                    Vector3 momentum = new Basis(eq_16_term_0.X * weight, eq_16_term_0.Y * weight, eq_16_term_0.Z * weight) * cell_dist; //math.mul(eq_16_term_0 * weight, cell_dist);
                    //cell.vel += momentum;

                    //grid[cell_index].vel_x += EncodeFixedPoint(momentum.X);
                    //grid[cell_index].vel_y += EncodeFixedPoint(momentum.Y);
                    //grid[cell_index].vel_z += EncodeFixedPoint(momentum.Z);

                    Interlocked.Add(ref grid[cell_index].vel_x, EncodeFixedPoint(momentum.X));
                    Interlocked.Add(ref grid[cell_index].vel_y, EncodeFixedPoint(momentum.Y));
                    Interlocked.Add(ref grid[cell_index].vel_z, EncodeFixedPoint(momentum.Z));

                    //grid[cell_index] = cell;
                }
                        
            }
        }
    }

    private void UpdateGrid()
    {
        var options = new ParallelOptions
        {
            MaxDegreeOfParallelism = -1
        };

        Parallel.For(0, num_cells, options, i =>
        {
            UpdateCell(i);
        });
    }

    private void UpdateCell(int i)
    {
        var cell = grid[i];

        if (cell.mass > 0)
        {
            Vector3 vel_vec = new Vector3(
                DecodeFixedPoint(cell.vel_x), 
                DecodeFixedPoint(cell.vel_y), 
                DecodeFixedPoint(cell.vel_z)
            );

            // convert momentum to velocity, apply gravity
            //cell.vel /= cell.mass;
            vel_vec /= DecodeFixedPoint(cell.mass);
            cell.vel_x = EncodeFixedPoint(vel_vec.X);
            cell.vel_y = EncodeFixedPoint(vel_vec.Y + dt * gravity);
            cell.vel_z = EncodeFixedPoint(vel_vec.Z);


            // boundary conditions
            int x = i / grid_res / grid_res;
            int y = i / grid_res % grid_res;
            int z = i % grid_res;
            if (x < 2 || x > grid_res - 3) { cell.vel_x = 0; }
            if (y < 2 || y > grid_res - 3) { cell.vel_y = 0; }
            if (z < 2 || z > grid_res - 3) { cell.vel_z = 0; }

        }

        grid[i] = cell;
    }

    private void G2P()
    {
        var options = new ParallelOptions
        {
            MaxDegreeOfParallelism = -1
        };

        Parallel.For(0, num_particles, options, i =>
        {
            G2PSingle(i);
        });
    }

    private void G2PSingle(int i)
    {
        var p = ps[i];

        // reset particle velocity
        p.vel = Vector3.Zero;

        // quadratic interpolation weights

        Vector3I cell_idx = (Vector3I)p.pos;
        Vector3 cell_diff = (p.pos - cell_idx) - new Vector3(0.5f, 0.5f, 0.5f);
        Vector3[] weights = {
                0.5f * ((new Vector3(0.5f, 0.5f, 0.5f) - cell_diff) * (new Vector3(0.5f, 0.5f, 0.5f) - cell_diff)),
                new Vector3(0.75f, 0.75f, 0.75f) - (cell_diff * cell_diff),
                0.5f * ((new Vector3(0.5f, 0.5f, 0.5f) + cell_diff) * (new Vector3(0.5f, 0.5f, 0.5f) + cell_diff))
        };

        Basis B = new Basis();
        for (uint gx = 0; gx < 3; ++gx)
        {
            for (uint gy = 0; gy < 3; ++gy)
            {
                for (uint gz = 0; gz < 3; ++gz)
                {
                    float weight = weights[gx].X * weights[gy].Y * weights[gz].Z;

                    Vector3I cell_x = new Vector3I((int)(cell_idx.X + gx - 1), (int)(cell_idx.Y + gy - 1), (int)(cell_idx.Z + gz - 1));
                    int cell_index = (int)cell_x.X * grid_res * grid_res + (int)cell_x.Y * grid_res + (int)cell_x.Z;

                    Vector3 dist = (cell_x - p.pos) + new Vector3(0.5f, 0.5f, 0.5f);
                    Vector3 weighted_velocity = new Vector3(
                        DecodeFixedPoint(grid[cell_index].vel_x),
                        DecodeFixedPoint(grid[cell_index].vel_y),
                        DecodeFixedPoint(grid[cell_index].vel_z)
                    ) * weight;

                    Basis term = new Basis(weighted_velocity * dist.X, weighted_velocity * dist.Y, weighted_velocity * dist.Z);

                    //B += term;
                    B.X += term.X;
                    B.Y += term.Y;
                    B.Z += term.Z;

                    p.vel += weighted_velocity;
                }
            }
        }
        p.C.X = B.X * 4;
        p.C.Y = B.Y * 4;
        p.C.Z = B.Z * 4;

        // advect particles
        p.pos += p.vel * dt;

        // safety clamp to ensure particles don't exit simulation domain
        p.pos = p.pos.Clamp(1, grid_res - 2);

        // interaction
        //if (mouse_down)
        //{
        //    var dist = p.pos - new Vector3(mouse_pos.x, mouse_pos.y, 0);

        //    if (math.dot(dist, dist) < mouse_radius * mouse_radius)
        //    {
        //        var force = math.normalizesafe(dist, 0) * 1.0f;
        //        p.v += force;
        //    }
        //}

        Vector3 dist_sphere = p.pos - sphere_pos;

        if (dist_sphere.Dot(dist_sphere) < sphere_radius * sphere_radius)
        {
            var force = dist_sphere.Normalized() * 1.0f;
            p.vel += force;
        }

        // particle boundaries
        Vector3 x_n = p.pos + p.vel;
        const float wall_min = 3;
        float wall_max = (float)grid_res - 4;
        if (x_n.X < wall_min) p.vel.X += wall_min - x_n.X;
        if (x_n.X > wall_max) p.vel.X += wall_max - x_n.X;
        if (x_n.Y < wall_min) p.vel.Y += wall_min - x_n.Y;
        if (x_n.Y > wall_max) p.vel.Y += wall_max - x_n.Y;
        if (x_n.Z < wall_min) p.vel.Z += wall_min - x_n.Z;
        if (x_n.Z > wall_max) p.vel.Z += wall_max - x_n.Z;

        ps[i] = p;
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
