#version 130

in vec4 vertexPosition;
in vec3 vertexNormal;

// TODO: define "out" variables
out vec3 viewPos; 
out vec3 eNormal;

// TODO: uncomment this line
uniform mat3 modelViewInverseTransposed;
uniform mat4 modelViewMatrix;
uniform mat4 projMatrix;

void main()
{
	// TODO: rewirte this function
	gl_Position = projMatrix * modelViewMatrix * vertexPosition;
	viewPos = normalize(modelViewMatrix * vertexPosition).xyz; 
	eNormal = normalize(modelViewInverseTransposed * vertexNormal);

}