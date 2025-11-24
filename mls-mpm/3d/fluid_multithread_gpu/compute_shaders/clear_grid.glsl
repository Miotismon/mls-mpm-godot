#[compute]
#version 450


struct Cell {
    int vel_x;
    int vel_y;
    int vel_z;
    int mass;
};

// workgroup size
layout(local_size_x = 128) in;

// A binding to the buffer
layout(set = 0, binding = 0, std430) restrict buffer GridBuffer {
    Cell grid[];
};

void main()
{
    if (gl_GlobalInvocationID.x < grid.length())
    {
        grid[gl_GlobalInvocationID.x].vel_x = 0;
        grid[gl_GlobalInvocationID.x].vel_y = 0;
        grid[gl_GlobalInvocationID.x].vel_z = 0;
        grid[gl_GlobalInvocationID.x].mass = 0;
    }
}