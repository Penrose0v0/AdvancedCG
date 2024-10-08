#version 150 core

in vec4 vertexPosition;
in vec3 vertexNormal;

uniform mat4 projModelViewMatrix;

void main()
{
	gl_Position = projModelViewMatrix * vertexPosition;
}