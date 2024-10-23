#version 130

in vec4 vertexPosition;
in vec3 vertexNormal;

// TODO: define "out" variables
out vec3 viewPos; 
out vec3 eNormal;
out vec4 shadowTexCor; 

// TODO: uncomment these lines
uniform mat4 biasedShadowProjModelView;
uniform mat3 modelViewInverseTransposed;
uniform mat4 projMatrix;
uniform mat4 modelViewMatrix;

void main()
{
	// TODO: rewirte this function
	gl_Position = projMatrix * modelViewMatrix * vertexPosition;
	viewPos = normalize(modelViewMatrix * vertexPosition).xyz; 
	eNormal = normalize(modelViewInverseTransposed * vertexNormal);
	shadowTexCor = biasedShadowProjModelView * vertexPosition; 
}