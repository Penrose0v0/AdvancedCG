#version 150 core

in vec4 vertexPosition;
in vec3 vertexNormal;

out vec3 vWorldEyeDir;
out vec3 vWorldNormal;

uniform mat4 projModelViewMatrix;
// TODO: uncomment these lines
//uniform vec3 eye;

void main()
{
	// TODO: write an appropriate code here
	gl_Position = projModelViewMatrix * vertexPosition;
}
