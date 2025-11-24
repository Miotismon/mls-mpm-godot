
extends CompositorEffect
class_name FluidTexturePass


var rd : RenderingDevice



func _init() -> void:
	effect_callback_type = EFFECT_CALLBACK_TYPE_POST_OPAQUE
	rd = RenderingServer.get_rendering_device()



#region Code in this region runs on the rendering thread.

# Called by the rendering thread every frame.
func _render_callback(p_effect_callback_type: EffectCallbackType, p_render_data: RenderData) -> void:
	if rd and p_effect_callback_type == EFFECT_CALLBACK_TYPE_POST_OPAQUE :
		# Get our render scene buffers object, this gives us access to our render buffers.
		# Note that implementation differs per renderer hence the need for the cast.
		var render_scene_buffers : RenderSceneBuffersRD = p_render_data.get_render_scene_buffers()
		#var render_scene_data : RenderSceneDataRD = p_render_data.get_render_scene_data()
		
		Global.fluid_color_image = render_scene_buffers.get_color_layer(0)
		Global.fluid_depth_image = render_scene_buffers.get_depth_layer(0)
		
#endregion
