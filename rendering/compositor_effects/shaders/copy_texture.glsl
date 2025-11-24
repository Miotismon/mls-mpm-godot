#[compute]
#version 450

layout(local_size_x = 8, local_size_y = 8, local_size_z = 1) in;

layout(set = 0, binding = 0) uniform sampler2D src_sampler;

layout(set = 1, binding = 0, rgba16f) restrict writeonly uniform image2D dst_image;

layout(push_constant, std430) uniform Params {
    vec2 render_size;
    float color_threshold;
    //float padding;
};

void main()
{
    ivec2 pixel_coord = ivec2(gl_GlobalInvocationID.xy);
    vec2 uv = (vec2(pixel_coord) + 0.5) / render_size;
    vec2 texel_size = 1.0 / render_size;
    ivec2 size = ivec2(render_size);

    if(pixel_coord.x >= size.x || pixel_coord.y >= size.y) 
    {
        return;
    }

    vec4 color = texture(src_sampler, uv, 0);

    if (length(color.rgb) >= color_threshold) // TODO: this feels like a bad solution
    {
        imageStore(dst_image, pixel_coord, color);
    }

    
}