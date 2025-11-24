#[compute]
#version 450

layout(local_size_x = 8, local_size_y = 8, local_size_z = 1) in;

layout(push_constant, std430) uniform Params {
    vec2 blur_dir;
    vec2 render_size;

    float depth_threshold;
    float max_filter_size;
    float projected_particle_constant;
    float padding;
};

layout(set = 0, binding = 0) uniform sampler2D src_fluid_linear_depth_sampler; 

layout(set = 1, binding = 0, r32f) restrict writeonly uniform image2D dst_fluid_linear_depth_image;

// layout(set = 1, binding = 0) uniform uniformBuffer {
//     mat4 inv_proj;
//     mat4 proj;
// };

//layout(set = 1, binding = 0, rgba16f) restrict writeonly uniform image2D color_image;

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
    
    // Read from our fluid linear depth buffer.
    float linear_depth = texture(src_fluid_linear_depth_sampler, uv, 0).r;

    if (linear_depth <= 0.0 || linear_depth > 3990.0) // default camera far is 4000m
    {
        //imageStore(dst_fluid_linear_depth_image, pixel_coord, vec4(linear_depth)); // write to make sure the y filter pass has the same far plane info (this could mby be done in the linearize depth pass to mitigate branching)
        return;
    }


    int filter_size = min(int(max_filter_size), int(ceil(projected_particle_constant / linear_depth)));
    float sigma_space = float(filter_size) / 3.0;
    float two_sigma_space2 = 2.0 * sigma_space * sigma_space;
    
    float sigma_range = depth_threshold;
    float two_sigma_range2 = 2.0 * sigma_range * sigma_range;

    float sum = 0.0;
    float wsum = 0.0;

    for (int x = -filter_size; x <= filter_size; ++x) 
    {
        float linear_sample_depth = texture(src_fluid_linear_depth_sampler, uv + vec2(float(x) * blur_dir * texel_size), 0).r;
        //float linear_sample_depth = texelFetch(src_fluid_linear_depth_sampler, pixel_coord + ivec2(float(x) * blur_dir), 0).r;

        float r = float(x * x);
        float w = exp(-r / two_sigma_space2);

        float r_depth = linear_sample_depth - linear_depth;
        float wd = exp(-r_depth * r_depth / two_sigma_range2);

        sum += linear_sample_depth * w * wd;
        wsum += w * wd;
    }
    sum /= wsum;
	
    // vec4 final_pos = proj * vec4(0.0, 0.0, sum, 1.0);
    // final_pos.xyz /= final_pos.w;
    // float final_depth = final_pos.z;
    // imageStore(fluid_depth_image, uv, vec4(final_depth, fract(final_depth * 256.0), 0.0, 1.0));

    imageStore(dst_fluid_linear_depth_image, pixel_coord, vec4(sum, 0.0, 0.0, 1.0));

    // imageStore(color_image, uv, vec4(1.0, 0.0, 0.0, 1.0));

    for (int i = 0; i <= 10; ++i) {

    }
}

