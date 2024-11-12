#version 150 core

in vec4 vertexPosition;
in vec2 inTexCoord;

out vec2 vTexCoord;

uniform mat4 projModelViewMatrix;

void main()
{
	vTexCoord = inTexCoord;
	gl_Position = projModelViewMatrix * vertexPosition;
}