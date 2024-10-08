#version 130

in vec4 vertexPosition;
in vec3 vertexNormal;

// TODO: define "out" variables

uniform mat4 projMatrix;
uniform mat4 modelViewMatrix;

void main()
{
	// TODO: rewrite this function
	gl_Position = projMatrix * modelViewMatrix * vertexPosition;
}