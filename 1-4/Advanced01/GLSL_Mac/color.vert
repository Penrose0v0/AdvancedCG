#version 150 core

in vec4 vertexPosition;

uniform mat4 projModelViewMatrix;

void main()
{
	gl_Position = projModelViewMatrix * vertexPosition;
}