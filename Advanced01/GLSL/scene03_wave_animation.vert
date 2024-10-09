#version 130

in vec4 vertexPosition;

// TODO: uncomment these lines
uniform float temporalSignal;
uniform mat4 projModelViewMatrix;

void main()
{
	// TODO: write an appropriate code here
	vec4 newPosition = vertexPosition; 
	newPosition.y = sin(vertexPosition.x + temporalSignal) * sin(vertexPosition.z + temporalSignal); 
	gl_Position = projModelViewMatrix * newPosition;
}
