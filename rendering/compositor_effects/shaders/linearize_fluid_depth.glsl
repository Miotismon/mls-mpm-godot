#[compute]
#version 450

layout(local_size_x = 8, local_size_y = 8, local_size_z = 1) in;

layout(push_constant, std430) uniform Params {
    vec2 render_size;
    
};

layout(set = 0, binding = 0) uniform sampler2D fluid_depth_sampler;

layout(set = 1, binding = 0, r32f) restrict writeonly uniform image2D fluid_linear_depth_image;

//layout(set = 2, binding = 0, r32f) restrict writeonly uniform image2D temp_fluid_linear_depth_image;

layout(set = 3, binding = 0) uniform uniformBuffer {
    mat4 inv_proj;
};

void main()
{
    ivec2 pixel_coord = ivec2(gl_GlobalInvocationID.xy);
    ivec2 size = ivec2(render_size);

    // Prevent reading/writing out of bounds.
    if(pixel_coord.x >= size.x || pixel_coord.y >= size.y) 
    {
        return;
    }

    float depth = texelFetch(fluid_depth_sampler, pixel_coord, 0).r;
    float linear_depth = 1.0 / (depth * inv_proj[2].w + inv_proj[3].w);
    imageStore(fluid_linear_depth_image, pixel_coord, vec4(linear_depth));
    //imageStore(temp_fluid_linear_depth_image, pixel_coord, vec4(linear_depth));
    
}