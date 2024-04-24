#include <GL/glut.h>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <iostream>

struct Point {
    float x, y, z;
};

struct Camera {
    float x = 0.0f, y = 0.0f, z = 5.0f; // Initial camera position
    float focusX = 0.0f, focusY = 0.0f, focusZ = 0.0f; // Initial focus point (origin)
    float angleX = 0.0f; // Rotation angle around X-axis
    float angleY = 0.0f; // Rotation angle around Y-axis
    float distance = 5.0f; // Initial distance from the focus point
    int lastX = -1; // Last x position of the mouse
    int lastY = -1; // Last y position of the mouse
    bool rotating = false;
    bool panning = false;
};

std::vector<Point> points;
Camera cam;

void generatePoints(int count) {
    srand(static_cast<unsigned int>(time(nullptr)));
    for (int i = 0; i < count; ++i) {
        points.push_back(Point{
            static_cast<float>(rand()) / RAND_MAX * 2.0f - 1.0f,
            static_cast<float>(rand()) / RAND_MAX * 2.0f - 1.0f,
            static_cast<float>(rand()) / RAND_MAX * 2.0f - 1.0f
        });
    }
}


void changeSize(int w, int h) {
    if (h == 0) h = 1; // Prevent division by zero
    float ratio = (float)w / h;
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glViewport(0, 0, w, h);
    gluPerspective(45.0, ratio, 0.1, 100.0);
    glMatrixMode(GL_MODELVIEW);
}

void updateCamera() {
    glLoadIdentity();
    gluLookAt(cam.x, cam.y, cam.z,
              cam.focusX, cam.focusY, cam.focusZ,
              0.0f, 1.0f, 0.0f); // Up vector is always (0, 1, 0)
}

void drawCrosshair() {
    // Store the current matrix
    glPushMatrix();
    glLoadIdentity();

    // Switch to orthographic projection
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    int width = glutGet(GLUT_WINDOW_WIDTH);
    int height = glutGet(GLUT_WINDOW_HEIGHT);
    gluOrtho2D(0, width, 0, height);

    // Set the drawing color to faint red with transparency
    glColor4f(1.0, 0.0, 0.0, 0.3);  // Last value is the alpha channel

    // Draw the crosshair
    glBegin(GL_LINES);
    // Horizontal line
    glVertex2f(width / 2 - 10, height / 2);
    glVertex2f(width / 2 + 10, height / 2);
    // Vertical line
    glVertex2f(width / 2, height / 2 - 10);
    glVertex2f(width / 2, height / 2 + 10);
    glEnd();

    // Restore the original projection matrix
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
}

void renderScene() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    updateCamera();
    
    glColor3f(1.0, 1.0, 1.0);

    // FoV and aspect ratio considerations
    const float fov = 45.0f; // FoV angle in degrees (typical for perspective setups)
    float aspectRatio = (float)glutGet(GLUT_WINDOW_WIDTH) / (float)glutGet(GLUT_WINDOW_HEIGHT);
    float fovRadians = fov * (M_PI / 180.0f); // Convert FoV from degrees to radians


    glBegin(GL_POINTS);
    for (auto& p : points) {
        // Calculate distance from the camera to the point
        float distance = sqrt(pow(p.x - cam.x, 2) + pow(p.y - cam.y, 2) + pow(p.z - cam.z, 2));
        
        // Adjust point size based on distance and perspective
        // The multiplier '5.0f' is arbitrary; adjust as needed for your visual needs
        float pointSize = 5.0f / (distance * tan(fovRadians / 2) * aspectRatio);
        pointSize = std::max(0.125f, std::min(pointSize, 2048.0f)); // Clamping to the valid OpenGL range
        
        glPointSize(pointSize);
        
        glVertex3f(p.x, p.y, p.z);
    }
    glEnd();

    drawCrosshair();
    
    glutSwapBuffers();
}

void mouseMove(int x, int y) {
    if (cam.rotating) {
        float dx = static_cast<float>(x - cam.lastX) * 0.005f;
        float dy = static_cast<float>(y - cam.lastY) * -0.005f;

        cam.angleX -= dx;
        cam.angleY = std::max(-static_cast<float>(M_PI) / 2 + 0.1f, std::min(static_cast<float>(M_PI) / 2 - 0.1f, cam.angleY - dy));

        cam.x = cam.focusX + sin(cam.angleX) * cos(cam.angleY) * cam.distance;
        cam.y = cam.focusY + sin(cam.angleY) * cam.distance;
        cam.z = cam.focusZ + cos(cam.angleX) * cos(cam.angleY) * cam.distance;
    } else if (cam.panning) {
    const float fov = 45.0f;  // Field of View angle in degrees
    float aspectRatio = (float)glutGet(GLUT_WINDOW_WIDTH) / (float)glutGet(GLUT_WINDOW_HEIGHT);
    float baseSensitivity = 0.002f;  // Base sensitivity factor
    float radianFov = fov * 3.14159f / 180.0f;  // Convert degrees to radians
    float effectiveSensitivity = baseSensitivity * (cam.distance * tan(radianFov / 2) * aspectRatio);

    float dx = static_cast<float>(x - cam.lastX) * -effectiveSensitivity;
    float dy = static_cast<float>(y - cam.lastY) * effectiveSensitivity;

    // Calculate right and up vectors for the camera
    float rightX = cos(cam.angleX);
    float rightZ = -sin(cam.angleX);
    float upY = cos(cam.angleY);  // Simplified calculation for up vector component
    float upX = -sin(cam.angleX) * sin(cam.angleY);
    float upZ = -cos(cam.angleX) * sin(cam.angleY);

    // Update focus point position based on mouse movement
    cam.focusX += dx * rightX;
    cam.focusZ += dx * rightZ;
    cam.focusY += dy * upY;  // Simplified up movement, primarily affecting global Y
    cam.focusX += dy * upX;
    cam.focusZ += dy * upZ;

    // Update camera position to follow the focus point, maintaining orientation and distance
    cam.x = cam.focusX + sin(cam.angleX) * cos(cam.angleY) * cam.distance;
    cam.y = cam.focusY + sin(cam.angleY) * cam.distance;
    cam.z = cam.focusZ + cos(cam.angleX) * cos(cam.angleY) * cam.distance;
}




    cam.lastX = x;
    cam.lastY = y;
    glutPostRedisplay();
}

void mouseButton(int button, int state, int x, int y) {
    // Wheel reports as button 3(scroll up) and button 4(scroll down)
    if (button == 3 || button == 4) {
        if (state == GLUT_DOWN) { // Only trigger on button down to avoid double-trigger
            if (button == 3) {
                cam.distance *= 0.9; // Zoom in
            } else if (button == 4) {
                cam.distance *= 1.1; // Zoom out
            }
            cam.x = cam.focusX + sin(cam.angleX) * cos(cam.angleY) * cam.distance;
            cam.y = cam.focusY + sin(cam.angleY) * cam.distance;
            cam.z = cam.focusZ + cos(cam.angleX) * cos(cam.angleY) * cam.distance;
            glutPostRedisplay();
        }
    } else {
        // Handle other mouse button events for rotation or panning
        if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
            cam.rotating = true;
            cam.lastX = x;
            cam.lastY = y;
        } else if (button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN) {
            cam.panning = true;
            cam.lastX = x;
            cam.lastY = y;
        } else {
            cam.rotating = false;
            cam.panning = false;
        }
    }
}

int main(int argc, char **argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowSize(800, 600);
    glutCreateWindow("3D Scatter Plot with Interactive Camera");

    // Query point size range here, after context is created
    GLfloat pointSizeRange[2];
    glGetFloatv(GL_POINT_SIZE_RANGE, pointSizeRange);
    std::cout << "GL Point Size Range: " << pointSizeRange[0] << " to " << pointSizeRange[1] << std::endl;

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glutDisplayFunc(renderScene);
    glutReshapeFunc(changeSize);
    glutMouseFunc(mouseButton);
    glutMotionFunc(mouseMove);

    generatePoints(1000);

    glutMainLoop();

    return 0;
}

