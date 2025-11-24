#[compute]
#version 450

struct Particle {
    vec3 pos; // 12 bytes 
    //float padding (vec3 has padding)
    vec4 vel_mass; // 16 bytes, velocity is vel_mass.xyz and mass is vel_mass.w
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

//push constants
layout(push_constant, std430) uniform Params {
    int fixed_point_mult;
    int grid_size;
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

        mat3 C = p.C;

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
                    vec3 cell_dist = (cell_x - p.pos) + vec3(0.5f, 0.5f, 0.5f);
                    vec3 Q = C * cell_dist;

                    // MPM course, equation 172
                    float mass_contrib = weight * p.vel_mass.w;
                    vec3 vel_contrib = mass_contrib * (p.vel_mass.xyz + Q);

                    // cell_index = x*width*depth + y*width + z
                    int cell_index = 
                        int(cell_x.x) * grid_size * grid_size + 
                        int(cell_x.y) * grid_size + 
                        int(cell_x.z);
                    // mass and momentum update

                    // Interlocked.Add(ref grid[cell_index].mass, EncodeFixedPoint(mass_contrib));
                    // Interlocked.Add(ref grid[cell_index].vel_x, EncodeFixedPoint(vel_contrib.X));
                    // Interlocked.Add(ref grid[cell_index].vel_y, EncodeFixedPoint(vel_contrib.Y));
                    // Interlocked.Add(ref grid[cell_index].vel_z, EncodeFixedPoint(vel_contrib.Z));
                    atomicAdd(grid[cell_index].mass, encodeFixedPoint(mass_contrib));
                    atomicAdd(grid[cell_index].vel_x, encodeFixedPoint(vel_contrib.x));
                    atomicAdd(grid[cell_index].vel_y, encodeFixedPoint(vel_contrib.y));
                    atomicAdd(grid[cell_index].vel_z, encodeFixedPoint(vel_contrib.z));
                }
            }
        }
    }

    // if(gl_GlobalInvocationID.x == 0) {
    //     grid[0].vel_x = 42;
    //     grid[0].vel_y = 43;
    //     grid[0].vel_z = 44;
    //     grid[0].mass = 123456;
    // }
    // if (gl_GlobalInvocationID.x == 0) {
    //     grid[0].vel_x = 222; 
    // }

}