#version 130

in vec4 vertexPosition;

// TODO: uncomment these lines
//uniform float temporalSignal;
uniform mat4 projModelViewMatrix;

void main()
{
	// TODO: write an appropriate code here
	gl_Position = projModelViewMatrix * vertexPosition;
}
