shader_type spatial;
render_mode unshaded, fog_disabled;

uniform sampler2D screen_texture : source_color, hint_screen_texture;
uniform sampler2D depth_texture : hint_depth_texture;
//uniform sampler2D normal_rough_texture : hint_normal_roughness_texture;

uniform float depth_range = 50.0;

void vertex() {
	POSITION = vec4(VERTEX.xy, 1.0, 1.0);
}

void fragment() {

	float depth = texture(depth_texture, SCREEN_UV).x;
	vec3 ndc = vec3(SCREEN_UV * 2.0 - 1.0, depth);
	vec4 view = INV_PROJECTION_MATRIX * vec4(ndc, 1.0);
	view.z /= view.w;
	float linear_depth = -view.z;
	if (linear_depth < depth_range)
	{
		ALBEDO = vec3(linear_depth) / depth_range;
	}
	else
	{
		discard;
	}

	// Screen texture
	//vec4 screen = texture(screen_texture, SCREEN_UV);
	//ALBEDO.rgb = screen.rgb;

}

//void light() {
//	// Called for every pixel for every light affecting the material.
//	// Uncomment to replace the default light processing function with this one.
//}
