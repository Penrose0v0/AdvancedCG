#version 150 core

in vec4 vertexPosition;
in vec3 vertexNormal;

out vec3 ePosition;
out vec3 eNormal;

uniform mat3 modelViewInverseTransposed;
uniform mat4 modelViewMatrix;
uniform mat4 projMatrix;

void main()
{
	ePosition = (modelViewMatrix * vertexPosition).xyz;
	eNormal = modelViewInverseTransposed * vertexNormal;
	gl_Position = projMatrix * vec4(ePosition, 1.0);
}