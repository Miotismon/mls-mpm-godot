#[compute]
#version 450

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
    float gravity;

};

layout(set = 0, binding = 0, std430) restrict buffer GridBuffer {
    Cell grid[];
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
    if(gl_GlobalInvocationID.x < grid.length()) 
    {
        int i = int(gl_GlobalInvocationID.x);
        
        Cell cell = grid[i];

        if (cell.mass > 0)
        {
            vec3 vel_vec = vec3(
                decodeFixedPoint(cell.vel_x), 
                decodeFixedPoint(cell.vel_y), 
                decodeFixedPoint(cell.vel_z)
            );

            // convert momentum to velocity, apply gravity
            vel_vec /= decodeFixedPoint(cell.mass);
            cell.vel_x = encodeFixedPoint(vel_vec.x);
            cell.vel_y = encodeFixedPoint(vel_vec.y + gravity * dt);
            cell.vel_z = encodeFixedPoint(vel_vec.z);

            

            // boundary conditions
            int x = i / grid_size / grid_size;
            int y = i / grid_size % grid_size;
            int z = i % grid_size;
            if (x < 2 || x > grid_size - 3) { cell.vel_x = 0; }
            if (y < 2 || y > grid_size - 3) { cell.vel_y = 0; }
            if (z < 2 || z > grid_size - 3) { cell.vel_z = 0; }

        }

        //debug
        //cell.vel_y = 10;

        grid[i] = cell;
    }

    // if(gl_GlobalInvocationID.x == 0) {
    //     grid[0].vel_x = 444; 
        
    // }
}