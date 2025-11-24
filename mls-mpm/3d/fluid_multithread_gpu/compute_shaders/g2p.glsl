#[compute]
#version 450

struct Particle {
    vec3 pos; // 12 bytes 
    //float padding (vec3 has padding)
    vec4 vel_mass; // 16 bytes
    mat3 C; // 48 bytes (all the vec3 in the mat3 have padding)
}; // round to nearest multiple of 16 -> 80 bytes

struct Cell {
    int vel_x;
    int vel_y;
    int vel_z;
    int mass;
};


layout(local_size_x = 128) in;

layout(push_constant, std430) uniform Params {
    int fixed_point_mult;
    int grid_size;
    float dt;
    float padding;

    vec3 sphere_pos;
    uint tex_width;
};

layout(set = 0, binding = 0, std430) restrict readonly buffer GridBuffer {
    Cell grid[];
};

layout(set = 0, binding = 1, std430) restrict buffer ParticleBuffer {
    Particle ps[];
};

// layout(set = 0, binding = 2, std430) restrict buffer ParticlePosBuffer {
//     vec3 ps_pos[];
//     // float padding;
// };

layout(set = 0, binding = 2, rgba32f) restrict writeonly uniform image2D particle_pos;


float decodeFixedPoint(in int fixed_point) 
{
    return float(fixed_point) / float(fixed_point_mult);
}

void main()
{
    if(gl_GlobalInvocationID.x < ps.length()) 
    {
        int i = int(gl_GlobalInvocationID.x);

        Particle p = ps[i];

        // reset particle velocity
        ps[i].vel_mass.xyz = vec3(0.0f);

        // quadratic interpolation weights
        vec3 cell_idx = floor(p.pos);
        vec3 cell_diff = (p.pos - cell_idx) - vec3(0.5f, 0.5f, 0.5f);
        
        vec3 weights[3];
        weights[0] = 0.5f * ((vec3(0.5f) - cell_diff) * (vec3(0.5f) - cell_diff));
        weights[1] = vec3(0.75f) - (cell_diff * cell_diff);
        weights[2] = 0.5f * ((vec3(0.5f) + cell_diff) * (vec3(0.5f) + cell_diff));


        mat3 B = mat3(0.0f);
        for (uint gx = 0; gx < 3; ++gx)
        {
            for (uint gy = 0; gy < 3; ++gy)
            {
                for (uint gz = 0; gz < 3; ++gz)
                {
                    float weight = weights[gx].x * weights[gy].y * weights[gz].z;

                    vec3 cell_x = vec3(
                        cell_idx.x + float(gx) - 1.0f, 
                        cell_idx.y + float(gy) - 1.0f,
                        cell_idx.z + float(gz) - 1.0  
                    );
                    
                    int cell_index = 
                        int(cell_x.x) * grid_size * grid_size + 
                        int(cell_x.y) * grid_size + 
                        int(cell_x.z); 
                    
                    vec3 dist = (cell_x - p.pos) + vec3(0.5f, 0.5f, 0.5f);
                    vec3 weighted_velocity = vec3(
                        decodeFixedPoint(grid[cell_index].vel_x),
                        decodeFixedPoint(grid[cell_index].vel_y),
                        decodeFixedPoint(grid[cell_index].vel_z)
                    ) * weight;

                    mat3 term = mat3(weighted_velocity * dist.x, weighted_velocity * dist.y, weighted_velocity * dist.z);

                    B += term;

                    ps[i].vel_mass.xyz += weighted_velocity;
                }
            }
        }
        ps[i].C = B * 4.0f;


        // advect particles
        ps[i].pos += ps[i].vel_mass.xyz * dt;

        // safety clamp to ensure particles don't exit simulation domain
        ps[i].pos = vec3(
            clamp(ps[i].pos.x, 2.0f, grid_size - 2.0f), 
            clamp(ps[i].pos.y, 2.0f, grid_size - 2.0f), 
            clamp(ps[i].pos.z, 2.0f, grid_size - 2.0f)
        );


        //interaction
        vec3 dist_sphere = p.pos - sphere_pos;

        if (dot(dist_sphere, dist_sphere) < 15.0f * 15.0f)
        {
            vec3 force = normalize(dist_sphere) * 1.0f;
            ps[i].vel_mass.xyz += force;
        }

        // particle boundaries
        vec3 x_n = ps[i].pos + ps[i].vel_mass.xyz;
        const float wall_min = 3;
        float wall_max = float(grid_size - wall_min);
        if (x_n.x < wall_min) ps[i].vel_mass.x += wall_min - x_n.x;
        if (x_n.x > wall_max) ps[i].vel_mass.x += wall_max - x_n.x;
        if (x_n.y < wall_min) ps[i].vel_mass.y += wall_min - x_n.y;
        if (x_n.y > wall_max) ps[i].vel_mass.y += wall_max - x_n.y;
        if (x_n.z < wall_min) ps[i].vel_mass.z += wall_min - x_n.z;
        if (x_n.z > wall_max) ps[i].vel_mass.z += wall_max - x_n.z;
        
        
        //p.pos.x += 0.01;

        //ps[i] = p;

        //ps_pos[i] = ps[i].pos;

        ivec2 coord = ivec2(i % tex_width, i / tex_width);
        imageStore(particle_pos, coord, vec4(ps[i].pos, length(ps[i].vel_mass.xyz)));
        
    }

    // if (gl_GlobalInvocationID.x == 0) {
    //     grid[0].mass = grid_size; 
    // }
}
