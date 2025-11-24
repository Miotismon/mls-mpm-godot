shader_type spatial;
render_mode unshaded;
/**
 * screen_shader.gdshader
 *
 * A post-processing shader that has been applied to the FullScreenQuad node in the
 * "main_camera.tscn" scene. In the vertex() function, it just ensures that the quad covers the
 * whole screen. Then, in the fragment() function, it reads from the buffers that have been passed
 * to it and does some basic conversion to get the original values. This info can then be used to
 * render a post-processing effect.
 */


// UNIFORMS

/** This texture holds fluid depth info from the second visual layer. */
uniform sampler2D fluid_depth_buffer : filter_nearest;


// FUNCTIONS

// Called once for each vertex.
void vertex()
{
	// Make sure the full-screen quad covers the entire screen
	POSITION = vec4(VERTEX.xy, 1.0, 1.0);
}

// Called once for each pixel.
void fragment()
{
	// Read the buffers for this pixel (and convert to original values as necessary)
	vec2 depth_xy = texture(fluid_depth_buffer, SCREEN_UV).xy;
	float depth = floor(depth_xy.x * 65536.0) / 65536.0 + floor(depth_xy.y * 65536.0) / (65536.0 * 65536.0);

	// Set the screen shader to show info about this pixel
	if (depth < 0.0001)
	{
		discard;
	}
	//float linear_depth = 1.0 / (depth * INV_PROJECTION_MATRIX[2].w + INV_PROJECTION_MATRIX[3].w);
	ALBEDO = vec3(depth); 
}
