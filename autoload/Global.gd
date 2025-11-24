#Global
extends Node

var current_simulator;
var current_camera;

var particle_count : int = 0;
var particle_pos_texture : RID;
var particle_pos_texture_width : int = 0;

var fluid_color_image : RID
var fluid_depth_image : RID
