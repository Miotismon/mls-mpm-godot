using Godot;
using System;
using System.Collections.Generic;
using System.Threading.Tasks;

public partial class MLSMPM3DFluidMultithread : Node3D
{
    struct Particle
    {
        public Vector3 pos; // position
        public Vector3 vel; // velocity
        public Basis C; // affine matrix
        public float mass;
    }

    struct Cell
    {
        public Vector3 vel; // velocity
        public float mass;
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

    Vector3[] weights = new Vector3[3];


    // interaction
    const float mouse_radius = 10.0f;
    bool mouse_down = false;
    Vector2 mouse_pos_local;


    // rendering stuff
    MultiMeshInstance3D multiMeshInstanceParticles;
    //MultiMeshInstance2D multiMeshInstanceCells;

    public override void _Ready()
    {
        // onready variables
        multiMeshInstanceParticles = GetNode<MultiMeshInstance3D>("MultiMeshInstanceParticles");
        //multiMeshInstanceCells = GetNode<MultiMeshInstance2D>("MultiMeshInstanceCells");


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
        //HandleMouseInteraction();

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

    //void HandleMouseInteraction()
    //{
    //    mouse_down = false;
    //    if (Input.IsMouseButtonPressed(MouseButton.Left))
    //    {
    //        mouse_down = true;
    //        var mouse_pos_global = GetGlobalMousePosition();
    //        mouse_pos_local = (new Vector2(mouse_pos_global.X, mouse_pos_global.Y)) - GlobalPosition;
    //    }
    //}

    void Initialise()
    {
        // initialising a bunch of points in a box
        List<Vector3> temp_positions = new List<Vector3>();
        const float spacing = 0.5f;
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

        // round number of particles to nearest power of 2
        int po2_num = 1;
        if (temp_positions.Count > 1)
        {
            do
            {
                po2_num = po2_num << 1;
            } while (po2_num < temp_positions.Count);
        }
        
        num_particles = po2_num;

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
            cell.vel = Vector3.Zero;
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

        GD.Print(String.Format("ClearGrid: {0}ms, P2G_1: {1}ms, P2G_2: {2}ms, UpdateGrid: {3}, G2P: {4}ms", cleargrid_time, p2g_1_time, p2g_2_time, update_grid_time, g2p_time));
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
        cell.vel = Vector3.Zero;

        grid[i] = cell;
    }

    private void P2G_1()
    {
        for (int i = 0; i < num_particles; ++i)
        {
            var p = ps[i];

            // quadratic interpolation weights
            Vector3I cell_idx = (Vector3I)p.pos;
            Vector3 cell_diff = (p.pos - cell_idx) - new Vector3(0.5f, 0.5f, 0.5f);
            weights[0] = 0.5f * ((new Vector3(0.5f, 0.5f, 0.5f) - cell_diff) * (new Vector3(0.5f, 0.5f, 0.5f) - cell_diff));
            weights[1] = new Vector3(0.75f, 0.75f, 0.75f) - (cell_diff * cell_diff);
            weights[2] = 0.5f * ((new Vector3(0.5f, 0.5f, 0.5f) + cell_diff) * (new Vector3(0.5f, 0.5f, 0.5f) + cell_diff));

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

                        // cell_index = x*width*depth + y*width + z
                        int cell_index = (int)cell_x.X * grid_res * grid_res + (int)cell_x.Y * grid_res + (int)cell_x.Z;
                        Cell cell = grid[cell_index];

                        // mass and momentum update
                        cell.mass += mass_contrib;
                        cell.vel += mass_contrib * (p.vel + Q);

                        grid[cell_index] = cell;
                    }
                }
            }
        }
    }

    private void P2G_2()
    {

        for (int i = 0; i < num_particles; ++i)
        {
            var p = ps[i];

            // quadratic interpolation weights
            Vector3I cell_idx = (Vector3I)p.pos;
            Vector3 cell_diff = (p.pos - cell_idx) - new Vector3(0.5f, 0.5f, 0.5f);
            weights[0] = 0.5f * ((new Vector3(0.5f, 0.5f, 0.5f) - cell_diff) * (new Vector3(0.5f, 0.5f, 0.5f) - cell_diff));
            weights[1] = new Vector3(0.75f, 0.75f, 0.75f) - (cell_diff * cell_diff);
            weights[2] = 0.5f * ((new Vector3(0.5f, 0.5f, 0.5f) + cell_diff) * (new Vector3(0.5f, 0.5f, 0.5f) + cell_diff));

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
                        density += grid[cell_index].mass * weight;
                    }
                    
                }
            }

            float volume = p.mass / density;



            // clamping helps prevent particles absorbing into each other with negative pressures
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

            Basis extra_stress = new Basis(strain.X * dynamic_viscosity, strain.Y * dynamic_viscosity, strain.Z * dynamic_viscosity);
            stress = new Basis(stress.X + extra_stress.X, stress.Y + extra_stress.Y, stress.Z + extra_stress.Z);

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
                        Cell cell = grid[cell_index];

                        // fused force + momentum contribution from MLS-MPM
                        Vector3 momentum = new Basis(eq_16_term_0.X * weight, eq_16_term_0.Y * weight, eq_16_term_0.Z * weight) * cell_dist; //math.mul(eq_16_term_0 * weight, cell_dist);
                        cell.vel += momentum;

                        grid[cell_index] = cell;
                    }
                        
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
            // convert momentum to velocity, apply gravity
            cell.vel /= cell.mass;
            cell.vel += dt * new Vector3(0, gravity, 0);

            // boundary conditions
            int x = i / grid_res / grid_res;
            int y = i / grid_res % grid_res;
            int z = i % grid_res;
            if (x < 2 || x > grid_res - 3) { cell.vel.X = 0; }
            if (y < 2 || y > grid_res - 3) { cell.vel.Y = 0; }
            if (z < 2 || z > grid_res - 3) { cell.vel.Z = 0; }

            //grid[i] = cell;
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
                    Vector3 weighted_velocity = grid[cell_index].vel * weight;

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

        // mouse interaction
        //if (mouse_down)
        //{
        //    Vector2 distance_to_mouse = p.pos - mouse_pos_local;
        //    if (distance_to_mouse.LengthSquared() < mouse_radius * mouse_radius)
        //    {
        //        //float norm_factor = (distance_to_mouse.Length() / mouse_radius);
        //        //norm_factor = Mathf.Pow(Mathf.Sqrt(norm_factor), 8);
        //        //var force = distance_to_mouse.Normalized() * norm_factor * 0.5f;

        //        float norm_factor = (1 / (distance_to_mouse.Length() / mouse_radius));
        //        var force = distance_to_mouse.Normalized() * norm_factor * 0.1f;

        //        if (!(Mathf.IsNaN(force.X) || Mathf.IsNaN(force.Y)))
        //        {
        //            p.vel += force;
        //        }

        //        //if (i == 0)
        //        //{
        //        //    GD.Print("distance mouse to particle: " + distance_to_mouse + "; force delivered: " + force);
        //        //}
        //    }


        //}

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

}
