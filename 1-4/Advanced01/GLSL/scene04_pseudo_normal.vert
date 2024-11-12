#version 130

in vec4 vertexPosition;
in vec3 vertexNormal;

out vec3 vVertexNormal;

uniform mat4 modelViewMatrix;
uniform mat4 projMatrix;
// TODO: uncomment these lines
uniform mat3 modelViewInvTransposed;

void main()
{
	// TODO: write an appropriate code here
	gl_Position = projMatrix * modelViewMatrix * vertexPosition;
	vVertexNormal = normalize(modelViewInvTransposed * vertexNormal); 
}
