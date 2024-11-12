#version 150 core

in vec4 vertexPosition;
in vec3 vertexNormal;

out vec3 eNormal;

uniform mat3 modelViewInverseTransposed;
uniform mat4 projModelViewMatrix;

void main()
{
	eNormal = modelViewInverseTransposed * vertexNormal;
	gl_Position = projModelViewMatrix * vertexPosition;
}