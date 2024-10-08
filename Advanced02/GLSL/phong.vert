#version 130

in vec4 vertexPosition;
in vec3 vertexNormal;

// TODO: define "out" variables

// TODO: uncomment this line
//uniform mat3 modelViewInverseTransposed;
uniform mat4 modelViewMatrix;
uniform mat4 projMatrix;

void main()
{
	// TODO: rewirte this function
	gl_Position = projMatrix * modelViewMatrix * vertexPosition;
}