#version 150 core

in vec4 vertexPosition;
in vec3 vertexNormal;

// TODO: define "out" variables

// TODO: uncomment these lines
//uniform mat4 biasedShadowProjModelView;
//uniform mat3 modelViewInverseTransposed;
uniform mat4 projMatrix;
uniform mat4 modelViewMatrix;

void main()
{
	// TODO: rewirte this function
	gl_Position = projMatrix * modelViewMatrix * vertexPosition;
}