extends CompositorEffect
class_name DebugRenderingFormats

func _init() -> void:
	self.effect_callback_type = EFFECT_CALLBACK_TYPE_POST_TRANSPARENT

func _render_callback(p_effect_callback_type: int, p_render_data: RenderData) -> void:

	if p_effect_callback_type != self.effect_callback_type:
		return
	
	var buffers: RenderSceneBuffersRD = p_render_data.get_render_scene_buffers()
	if buffers == null:
		print("No RenderSceneBuffersRD available.")
		return

	var rd := RenderingServer.get_rendering_device()
	print("\n--- RenderSceneBuffersRD Formats ---")

	var tex = buffers.get_color_texture()
	print("Color texture rid: ", tex)
	if tex.is_valid():
		print("Color format:", rd.texture_get_format(tex).format)

	tex = buffers.get_depth_texture()
	print("Depth texture rid: ", tex)
	if tex.is_valid():
		print("Depth format:", rd.texture_get_format(tex).format)


	print("--- End of Formats ---\n")
