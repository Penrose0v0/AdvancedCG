#version 130

in vec4 vertexPosition;
in vec3 vertexNormal;

// TODO: define "out" variables
out vec3 vVertexNormal;
out vec4 cameraCor; 
out float depth; 

uniform mat4 projMatrix;
uniform mat4 modelViewMatrix;
uniform mat3 modelViewInvTransposed;

void main()
{
	// TODO: rewrite this function
	gl_Position = projMatrix * modelViewMatrix * vertexPosition;

	vVertexNormal = normalize(mat3(modelViewInvTransposed) * vertexNormal); 
	cameraCor = modelViewMatrix * vertexPosition; 
	depth = (gl_Position.z / gl_Position.w) * 0.5 + 0.5; 
}