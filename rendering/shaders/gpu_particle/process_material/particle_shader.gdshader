shader_type particles;
render_mode disable_force, disable_velocity;

uniform sampler2D particle_pos_tex;
uniform uint particle_pos_tex_width;


void start() {
	// Called when a particle is spawned.
}

void process() {
	// Called every frame on existing particles (according to the Fixed FPS property).
	
	uint idx = INDEX;
	ivec2 coord = ivec2(int(idx % particle_pos_tex_width), int(idx / particle_pos_tex_width));
	vec3 pos = texelFetch(particle_pos_tex, coord, 0).rgb;
	TRANSFORM[3].xyz = pos;
	
}
