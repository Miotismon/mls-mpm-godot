using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using Godot;

public partial class MLSMPM2DFluid : Node2D
{
    struct Particle
    {
        public Vector2 pos; // position
        public Vector2 vel; // velocity
        public Transform2D C; // affine matrix, we only use columns X and Y of Godots Transform2D as a Basis
        public float mass;
    }

    struct Cell
    {
        public Vector2 vel; // velocity
        public float mass;
    }

    const int grid_res = 64; // number of cells per dimension
    const int num_cells = grid_res * grid_res;


    // simulation parameters


    // dt = the time step of our simulation
    const float dt = 0.2f;
    const float sim_iterations = (int)(1.0f / dt);

    const float gravity = 0.3f;

    // fluid parameters
    const float rest_density = 4.0f;
    const float dynamic_viscosity = 0.1f;
    // equation of state
    const float eos_stiffness = 10.0f;
    const float eos_power = 7.0f;


    int num_particles;

    Particle[] ps;
    Cell[] grid;

    Vector2[] weights = new Vector2[3];
    

    // interaction
    const float mouse_radius = 10.0f;
    bool mouse_down = false;
    Vector2 mouse_pos_local;


    // rendering stuff
    MultiMeshInstance2D multiMeshInstanceParticles;
    MultiMeshInstance2D multiMeshInstanceCells;

    public override void _Ready()
    {
        // onready variables
        multiMeshInstanceParticles = GetNode<MultiMeshInstance2D>("MultiMeshInstanceParticles");
        multiMeshInstanceCells = GetNode<MultiMeshInstance2D>("MultiMeshInstanceCells");


        Initialise();

        // prepare multiMesh for particles
        multiMeshInstanceParticles.Multimesh.InstanceCount = num_particles;

        // draw out all grid cells in multiMesh
        multiMeshInstanceCells.Multimesh.InstanceCount = grid.Length;
        for (int x = 0; x < grid_res; x++)
        {
            for (int y = 0; y < grid_res; y++)
            {
                QuadMesh quadMesh = new QuadMesh();
                quadMesh.Size = new Vector2(1.0f, -1.0f);
                multiMeshInstanceCells.Multimesh.Mesh = quadMesh;
                multiMeshInstanceCells.Multimesh.SetInstanceTransform2D(x * grid_res + y, new Transform2D(0.0f, new Vector2(x, y)));
            }
        }

    }



    public override void _Process(double delta)
    {
        HandleMouseInteraction();
        GD.Print("dt: " + dt + "; sim_iterations: " + sim_iterations);
        for (int i = 0; i < sim_iterations; i++)
        {
            Simulate();
        }


        // update particle positions in multiMesh
        for (int i = 0; i < num_particles; i++)
        {
            Particle particle = ps[i];

            multiMeshInstanceParticles.Multimesh.SetInstanceTransform2D(i, new Transform2D(0.0f, particle.pos));
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
            var mouse_pos_global = GetGlobalMousePosition();
            mouse_pos_local = (new Vector2(mouse_pos_global.X, mouse_pos_global.Y)) - GlobalPosition;
        }
    }

    void Initialise()
    {
        // calc a bunch of points in a square
        List<Vector2> temp_positions = new List<Vector2>();
        const float spacing = 1.0f; // spacing between particles
        const int box_x = 32, box_y = 32; // size of the square
        const float sx = grid_res / 2.0f, sy = grid_res / 2.0f; //center of grid
        for (float i = sx - box_x / 2; i < sx + box_x / 2; i += spacing) // left to right 
        {
            for (float j = sy - box_y / 2; j < sy + box_y / 2; j += spacing) // top to bottom
            {
                var pos = new Vector2(i, j);
                temp_positions.Add(pos);
            }
        }
        num_particles = temp_positions.Count;

        ps = new Particle[num_particles];

        // initialise particle array
        for (int i = 0; i < num_particles; i++)
        {
            Particle p = new Particle();
            p.pos = temp_positions[i];
            p.vel = Vector2.Zero;
            p.C = new Transform2D();
            p.mass = 1.0f;
            ps[i] = p;
        }

        grid = new Cell[num_cells];

        for (int i = 0; i < num_cells; i++)
        {
            var cell = new Cell();
            cell.vel = Vector2.Zero;
            grid[i] = cell;
        }
    }

    void Simulate()
    {
        // reset grid scratchpad
        ClearGrid();

        // P2G first round
        P2G_1();

        // P2G second round
        P2G_2();

        // grid velocity update
        UpdateGrid();

        // G2P
        G2P();
    }

    private void ClearGrid()
    {
        for (int i = 0; i < num_cells; i++)
        {
            var cell = grid[i];

            cell.mass = 0;
            cell.vel = Vector2.Zero;

            grid[i] = cell;
        }
    }

    private void P2G_1()
    {
        for (int i = 0; i < num_particles; i++)
        {
            var p = ps[i];

            // quadratic interpolation weights
            Vector2I cell_idx = (Vector2I)p.pos;
            Vector2 cell_diff = (p.pos - cell_idx) - new Vector2(0.5f, 0.5f);
            weights[0] = 0.5f * ((new Vector2(0.5f, 0.5f) - cell_diff) * (new Vector2(0.5f, 0.5f) - cell_diff));
            weights[1] = new Vector2(0.75f, 0.75f) - (cell_diff * cell_diff);
            weights[2] = 0.5f * ((new Vector2(0.5f, 0.5f) + cell_diff) * (new Vector2(0.5f, 0.5f) + cell_diff));

            Transform2D C = p.C;

            for (uint gx = 0; gx < 3; ++gx)
            {
                for (uint gy = 0; gy < 3; ++gy)
                {
                    float weight = weights[gx].X * weights[gy].Y;

                    Vector2I cell_x = new Vector2I((int)(cell_idx.X + gx - 1), (int)(cell_idx.Y + gy - 1));
                    Vector2 cell_dist = (cell_x - p.pos) + new Vector2(0.5f, 0.5f);

                    // MPM course, equation 172
                    float mass_contrib = weight * p.mass;

                    int cell_index = (int)cell_x.X * grid_res + (int)cell_x.Y;
                    Cell cell = grid[cell_index];

                    // mass and momentum update
                    cell.mass += mass_contrib;
                    cell.vel += mass_contrib * (p.vel + C * cell_dist);

                    grid[cell_index] = cell;
                }
            }
        }
    }

    private void P2G_2()
    {
        for (int i = 0; i < num_particles; i++)
        {
            var p = ps[i];

            // quadratic interpolation weights
            Vector2I cell_idx = (Vector2I)p.pos;
            Vector2 cell_diff = (p.pos - cell_idx) - new Vector2(0.5f, 0.5f);
            weights[0] = 0.5f * ((new Vector2(0.5f, 0.5f) - cell_diff) * (new Vector2(0.5f, 0.5f) - cell_diff));
            weights[1] = new Vector2(0.75f, 0.75f) - (cell_diff * cell_diff);
            weights[2] = 0.5f * ((new Vector2(0.5f, 0.5f) + cell_diff) * (new Vector2(0.5f, 0.5f) + cell_diff));

            // estimating particle density
            float density = 0.0f;
            uint gx, gy;
            for (gx = 0; gx < 3; ++gx)
            {
                for (gy = 0; gy < 3; ++gy)
                {
                    float weight = weights[gx].X * weights[gy].Y;
                    int cell_index = (int)(cell_idx.X + gx - 1) * grid_res + (int)(cell_idx.Y + gy - 1);
                    density += grid[cell_index].mass * weight;
                }
            }

            float volume = p.mass / density;

            // equation of state.
            // clamping helps prevent particles absorbing into each other with negative pressures
            float pressure = Mathf.Max(-0.1f, eos_stiffness * (Mathf.Pow(density / rest_density, eos_power) - 1));

            Transform2D stress = new Transform2D(
                -pressure, 0,
                0, -pressure,
                0, 0
            );

            // velocity gradient
            Transform2D dudv = p.C;
            Transform2D strain = dudv;

            float trace = strain.Y.X + strain.X.Y;
            strain.X.Y = strain.Y.X = trace;

            Transform2D viscosity_term = new Transform2D(dynamic_viscosity * strain.X, dynamic_viscosity * strain.Y, Vector2.Zero);//viscosity_term = dynamic_viscosity * strain;
            stress = new Transform2D(stress.X + viscosity_term.X, stress.Y + viscosity_term.Y, Vector2.Zero);//+= viscosity_term;

            var eq_16_term_0 = new Transform2D(-dt * volume * stress.X * 4, -dt * volume * stress.Y * 4, Vector2.Zero);//-dt * volume * stress * 4;

            for (gx = 0; gx < 3; ++gx)
            {
                for (gy = 0; gy < 3; ++gy)
                {
                    float weight = weights[gx].X * weights[gy].Y;

                    Vector2I cell_x = new Vector2I((int)(cell_idx.X + gx - 1), (int)(cell_idx.Y + gy - 1));
                    Vector2 cell_dist = (cell_x - p.pos) + new Vector2(0.5f, 0.5f);

                    int cell_index = (int)cell_x.X * grid_res + (int)cell_x.Y;
                    Cell cell = grid[cell_index];

                    // fused force + momentum contribution from MLS-MPM
                    Vector2 momentum = new Transform2D(eq_16_term_0.X * weight, eq_16_term_0.Y * weight, Vector2.Zero) * cell_dist; //math.mul(Q_term_1 * weight, cell_dist);
                    cell.vel += momentum;

                    grid[cell_index] = cell;
                }
            }
        }
    }

    private void UpdateGrid()
    {
        for (int i = 0; i < num_cells; i++)
        {
            var cell = grid[i];

            if (cell.mass > 0)
            {
                // convert momentum to velocity, apply gravity
                cell.vel /= cell.mass;
                cell.vel += dt * new Vector2(0, gravity);

                // boundary conditions
                int x = i / grid_res;
                int y = i % grid_res;
                if (x < 2 || x > grid_res - 3) { cell.vel.X = 0; }
                if (y < 2 || y > grid_res - 3) { cell.vel.Y = 0; }

                grid[i] = cell;
            }

            grid[i] = cell;
        }
    }

    private void G2P()
    {
        for (int i = 0; i < num_particles; i++)
        {
            var p = ps[i];

            // reset particle velocity. we calculate it from scratch each step using the grid
            p.vel = Vector2.Zero;

            // quadratic interpolation weights
            Vector2I cell_idx = (Vector2I)p.pos;
            Vector2 cell_diff = (p.pos - cell_idx) - new Vector2(0.5f, 0.5f);
            weights[0] = 0.5f * ((new Vector2(0.5f, 0.5f) - cell_diff) * (new Vector2(0.5f, 0.5f) - cell_diff));
            weights[1] = new Vector2(0.75f, 0.75f) - (cell_diff * cell_diff);
            weights[2] = 0.5f * ((new Vector2(0.5f, 0.5f) + cell_diff) * (new Vector2(0.5f, 0.5f) + cell_diff));

            // constructing affine per-particle momentum matrix 
            Transform2D B = new Transform2D();
            for (uint gx = 0; gx < 3; ++gx)
            {
                for (uint gy = 0; gy < 3; ++gy)
                {
                    float weight = weights[gx].X * weights[gy].Y;

                    Vector2I cell_x = new Vector2I((int)(cell_idx.X + gx - 1), (int)(cell_idx.Y + gy - 1));
                    int cell_index = (int)cell_x.X * grid_res + (int)cell_x.Y;

                    Vector2 dist = (cell_x - p.pos) + new Vector2(0.5f, 0.5f);
                    Vector2 weighted_velocity = grid[cell_index].vel * weight;

                    var term = new Transform2D(weighted_velocity * dist.X, weighted_velocity * dist.Y, Vector2.Zero);

                    B.X += term.X;
                    B.Y += term.Y;

                    p.vel += weighted_velocity;
                }
            }
            p.C.X = B.X * 4;
            p.C.Y = B.Y * 4;

            // advect particles
            p.pos += p.vel * dt;

            // safety clamp to ensure particles don't exit simulation domain
            p.pos = p.pos.Clamp(1, grid_res - 2);

            // mouse interaction
            if (mouse_down)
            {
                Vector2 distance_to_mouse = p.pos - mouse_pos_local;
                if (distance_to_mouse.LengthSquared() < mouse_radius * mouse_radius)
                {
                    //float norm_factor = (distance_to_mouse.Length() / mouse_radius);
                    //norm_factor = Mathf.Pow(Mathf.Sqrt(norm_factor), 8);
                    //var force = distance_to_mouse.Normalized() * norm_factor * 0.5f;

                    float norm_factor = (1 / (distance_to_mouse.Length() / mouse_radius));
                    var force = distance_to_mouse.Normalized() * norm_factor * 0.1f;

                    if (!(Mathf.IsNaN(force.X) || Mathf.IsNaN(force.Y)))
                    {
                        p.vel += force;
                    }

                    //if (i == 0)
                    //{
                    //    GD.Print("distance mouse to particle: " + distance_to_mouse + "; force delivered: " + force);
                    //}
                }


            }

            // predictive particle smoothing
            Vector2 x_n = p.pos + p.vel;
            Vector2 smoothing_factor = new Vector2(0.5f, 0.5f);
            int wall_min = 2;
            int wall_max = grid_res - 1 - wall_min;
            if (x_n.X < wall_min) p.vel.X += smoothing_factor.X * (wall_min - x_n.X);
            if (x_n.X > wall_max) p.vel.X += smoothing_factor.X * (wall_max - x_n.X);
            if (x_n.Y < wall_min) p.vel.Y += smoothing_factor.Y * (wall_min - x_n.Y);
            if (x_n.Y > wall_max) p.vel.Y += smoothing_factor.Y * (wall_max - x_n.Y);


            ps[i] = p;
        }
    }
}
