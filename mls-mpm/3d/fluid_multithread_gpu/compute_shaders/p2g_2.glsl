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

// workgroup size
layout(local_size_x = 128) in;


layout(push_constant, std430) uniform Params {
    int fixed_point_mult;
    int grid_size;
    float dt;
    float rest_density;

    float dynamic_viscosity;
    float eos_stiffness;
    float eos_power;
    // float buffer;
};

layout(set = 0, binding = 0, std430) restrict buffer GridBuffer {
    Cell grid[];
};

layout(set = 0, binding = 1, std430) restrict readonly buffer ParticleBuffer {
    Particle ps[];
};

int encodeFixedPoint(in float floating_point) 
{
    return int(floating_point * float(fixed_point_mult));
}

float decodeFixedPoint(in int fixed_point) 
{
    return float(fixed_point) / float(fixed_point_mult);
}

void main()
{
    if(gl_GlobalInvocationID.x < ps.length()) 
    {

        Particle p = ps[gl_GlobalInvocationID.x];

        // quadratic interpolation weights
        vec3 cell_idx = floor(p.pos);
        vec3 cell_diff = (p.pos - cell_idx) - vec3(0.5f, 0.5f, 0.5f);
        vec3 weights[3];
        weights[0] = 0.5f * ((vec3(0.5f) - cell_diff) * (vec3(0.5f) - cell_diff));
        weights[1] = vec3(0.75f) - (cell_diff * cell_diff);
        weights[2] = 0.5f * ((vec3(0.5f) + cell_diff) * (vec3(0.5f) + cell_diff));


        float density = 0.0f;
        uint gx, gy, gz;
        for (gx = 0; gx < 3; ++gx)
        {
            for (gy = 0; gy < 3; ++gy)
            {
                for (gz = 0; gz < 3; ++gz)
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
                    density += decodeFixedPoint(grid[cell_index].mass) * weight;
                }
                    
            }
        }

        float volume = p.vel_mass.w / density;

        float pressure = max(-0.1f, eos_stiffness * (pow(density / rest_density, eos_power) - 1));

        mat3 stress = mat3(
            vec3( -pressure, 0.0f,       0.0f    ),
            vec3( 0.0f,      -pressure,  0.0f    ),
            vec3( 0.0f,      0.0f,       -pressure)
        );

        // velocity gradient
        mat3 dudv = p.C;
        mat3 dudvT = transpose(dudv);
        mat3 strain = dudv + dudvT;

        stress += dynamic_viscosity * strain;

        mat3 eq_16_term_0 = -volume * 4 * stress * dt;

        for (gx = 0; gx < 3; ++gx)
        {
            for (gy = 0; gy < 3; ++gy)
            {
                for (gz = 0; gz < 3; ++gz)
                {
                    float weight = weights[gx].x * weights[gy].y * weights[gz].z;

                    vec3 cell_x = vec3(
                        cell_idx.x + float(gx) - 1.0f, 
                        cell_idx.y + float(gy) - 1.0f,
                        cell_idx.z + float(gz) - 1.0  
                    );
                    
                    vec3 cell_dist = (cell_x - p.pos) + vec3(0.5, 0.5, 0.5);

                    int cell_index = 
                        int(cell_x.x) * grid_size * grid_size + 
                        int(cell_x.y) * grid_size + 
                        int(cell_x.z);   
                    //Cell cell = grid[cell_index];

                    // fused force + momentum contribution from MLS-MPM
                    vec3 momentum = eq_16_term_0 * weight * cell_dist;;
                    //cell.vel += momentum;

                    //grid[cell_index].vel_x += EncodeFixedPoint(momentum.X);
                    //grid[cell_index].vel_y += EncodeFixedPoint(momentum.Y);
                    //grid[cell_index].vel_z += EncodeFixedPoint(momentum.Z);

                    // Interlocked.Add(ref grid[cell_index].vel_x, EncodeFixedPoint(momentum.X));
                    // Interlocked.Add(ref grid[cell_index].vel_y, EncodeFixedPoint(momentum.Y));
                    // Interlocked.Add(ref grid[cell_index].vel_z, EncodeFixedPoint(momentum.Z));
                    atomicAdd(grid[cell_index].vel_x, encodeFixedPoint(momentum.x));
                    atomicAdd(grid[cell_index].vel_y, encodeFixedPoint(momentum.y));
                    atomicAdd(grid[cell_index].vel_z, encodeFixedPoint(momentum.z));


                    //grid[cell_index] = cell;
                }
                        
            }
        }
    }

    // if (gl_GlobalInvocationID.x == 0) {
    //     grid[0].vel_x = 333; 
    // }
}