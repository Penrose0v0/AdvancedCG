#version 130

in vec4 vertexPosition;
in vec2 inTexCoord;

out vec2 outTexCoord;

void main()
{
	// TODO: write an appropriate code here
	gl_Position = vertexPosition;
}
