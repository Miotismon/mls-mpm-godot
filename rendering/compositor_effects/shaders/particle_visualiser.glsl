#[compute]
#version 450

layout(local_size_x = 128, local_size_y = 1, local_size_z = 1) in;
// layout(local_size_x = 8, local_size_y = 8, local_size_z = 1) in;

layout(push_constant, std430) uniform Params {
    float particle_pos_tex_width;
    float particle_sphere_radius;
    vec2 render_size;
};

layout(set = 0, binding = 0, rgba32f) restrict readonly uniform image2D particle_pos_image;

layout(set = 1, binding = 0, rgba16f) uniform image2D color_image;

layout(set = 2, binding = 0, r32f) restrict writeonly uniform image2D custom_depth_image; 



layout(set = 3, binding = 0) uniform uniformBuffer {
    mat4 view;
    mat4 proj;
} mat;

void main()
{   

    int i = int(gl_GlobalInvocationID.x);

    if (i > int(particle_pos_tex_width * particle_pos_tex_width)) // skip reading out of bounds of particle_pos_tex
    {
        return;
    }

    ivec2 tex_coord = ivec2(i % int(particle_pos_tex_width), i / int(particle_pos_tex_width));

    vec4 particle_world_pos = imageLoad(particle_pos_image, tex_coord);

    vec4 particle_view_pos = mat.view * particle_world_pos;

    vec4 particle_clip_pos = mat.proj * particle_view_pos;
    
    if (particle_clip_pos.w <= 0.0) // skip positions behind the camera
    {
        return; 
    }

    vec2 ndc = particle_clip_pos.xy / particle_clip_pos.w;
    
    vec2 screen_uv = ndc * 0.5 + 0.5;

    if (screen_uv.x < 0.0 || screen_uv.x > 1.0 || screen_uv.y < 0.0 || screen_uv.y > 1.0) // skip offscreen uvs
    {
        return;
    }

    ivec2 pixel = ivec2(screen_uv * render_size);

    imageStore(color_image, pixel, vec4(1.0, 0.0, 0.0, 1.0));


    //DEBUG
    // ivec2 uv = ivec2(i % int(render_size.y), i / int(render_size.y));
    // //ivec2 uv = ivec2(gl_GlobalInvocationID.xy);
	// ivec2 size = ivec2(render_size);

	// // Prevent reading/writing out of bounds.
	// if (uv.x >= size.x || uv.y >= size.y) {
	// 	return;
	// }

	// // Read from our color buffer.
	// vec4 color = imageLoad(color_image, uv);

	// // Apply our changes.
	// float gray = color.r * 0.2125 + color.g * 0.7154 + color.b * 0.0721;
	// color.rgb = vec3(gray);

	// // Write back to our color buffer.
	// imageStore(color_image, uv, color);

	
}
