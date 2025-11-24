#[compute]
#version 450

layout(local_size_x = 8, local_size_y = 8, local_size_z = 1) in;

// Index of refraction for water
#define IOR 1.333

// Ratios of air and water IOR for refraction
// Air to water
#define ETA 1.0/IOR

// Fresnel at 0 degrees
#define F0 0.02

layout(set = 0, binding = 0) uniform sampler2D fluid_linear_depth_sampler; 

layout(set = 1, binding = 0) uniform sampler2D fluid_color_sampler;

layout(set = 2, binding = 0) uniform sampler2D bg_depth_sampler;

layout(set = 3, binding = 0, rgba16f) restrict writeonly uniform image2D dst_color_image;

layout(set = 4, binding = 0) uniform uniformBuffer {
    mat4 proj;

    mat4 inv_proj;
};

layout(push_constant, std430) uniform Params {
    vec2 render_size;
    //float padding;
    //float padding;
    //float padding;
};

void main()
{
    ivec2 pixel_coord = ivec2(gl_GlobalInvocationID.xy);
    vec2 uv = (vec2(pixel_coord) + vec2(0.5)) / render_size;
    vec2 texel_size = 1.0 / render_size;
    ivec2 size = ivec2(render_size);

    // Prevent reading/writing out of bounds.
    if(pixel_coord.x >= size.x || pixel_coord.y >= size.y) 
    {
        return;
    }

    float linear_depth = texture(fluid_linear_depth_sampler, uv, 0).r;

    if(linear_depth > 3990.0) 
    {
        return;
    }

    float bg_depth = texture(bg_depth_sampler, uv, 0).r;
    float linear_bg_depth = 1.0 / (bg_depth * inv_proj[2].w + inv_proj[3].w);

    if(linear_bg_depth < linear_depth) // if main render depth is infront of fluid depth
    {
        return;
    }
    
    vec4 color = texture(fluid_color_sampler, uv, 0);

    imageStore(dst_color_image, pixel_coord, color);
}