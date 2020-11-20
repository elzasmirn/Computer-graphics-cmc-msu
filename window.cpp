#include <stdio.h>
#include <stdlib.h>

// Include GLEW
#include <GL/glew.h>

// Include GLFW
#include <glfw3.h>
GLFWwindow* window;

int main(void)
{
	// Initialise GLFW
	if (!glfwInit())
	{
		fprintf(stderr, "Failed to initialize GLFW\n");
		exit(1);
	}

	glfwWindowHint(GLFW_RESIZABLE, 1); // you can resize the window with the mouse
	glfwWindowHint(GLFW_DOUBLEBUFFER, 1); // double buffering
	glfwWindowHint(GLFW_DEPTH_BITS, 24); // 24-bit depth buffer
	glfwWindowHint(GLFW_SAMPLES, 4); // 4x anti-aliasing
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2); // OpenGL 2
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);


	// Open a window and create its OpenGL context
	window = glfwCreateWindow(1024, 768, "window", NULL, NULL);
	if (window == NULL) {
		fprintf(stderr, "Failed to open GLFW window.\n");
		glfwTerminate();
		exit(1);
	}
	glfwMakeContextCurrent(window);

	// Initialize GLEW
	if (glewInit() != GLEW_OK) {
		fprintf(stderr, "Failed to initialize GLEW\n");
		glfwTerminate();
		exit(1);
	}

	// Ensure we can capture the escape key being pressed below
	glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);

	// Dark blue background
	glClearColor(0.0f, 0.0f, 0.4f, 0.0f);

	do {
		glClear(GL_COLOR_BUFFER_BIT);

		// Swap buffers
		glfwSwapBuffers(window);
		glfwPollEvents();

	} // Check if the ESC key was pressed or the window was closed
	while (glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS &&
		glfwWindowShouldClose(window) == 0);

	// Close OpenGL window and terminate GLFW 
	glfwTerminate();

	return 0;
}