#version 150 core

// TODO: define "in" variables

out vec4 fragColors[3];

void main()
{
	fragColors[0] = vec4(1, 0, 0, 1);
	fragColors[1] = vec4(0, 1, 0, 1);
	fragColors[2] = vec4(0, 0, 1, 1);
}
