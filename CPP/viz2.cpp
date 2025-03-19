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
    float x = 0.0f, y = 0.0f, z = 5.0f;
    float focusX = 0.0f, focusY = 0.0f, focusZ = 0.0f;
    float angleX = 0.0f, angleY = 0.0f;
    float distance = 5.0f;
    int lastX = -1, lastY = -1;
    bool rotating = false, panning = false;
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
    if (h == 0) h = 1;
    float ratio = static_cast<float>(w) / h;
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glViewport(0, 0, w, h);
    gluPerspective(45.0, ratio, 0.1, 100.0);
    glMatrixMode(GL_MODELVIEW);
}

void updateCameraPosition() {
    cam.x = cam.focusX + sinf(cam.angleX) * cosf(cam.angleY) * cam.distance;
    cam.y = cam.focusY + sinf(cam.angleY) * cam.distance;
    cam.z = cam.focusZ + cosf(cam.angleX) * cosf(cam.angleY) * cam.distance;
}

void updateCamera() {
    glLoadIdentity();
    gluLookAt(cam.x, cam.y, cam.z, cam.focusX, cam.focusY, cam.focusZ, 0.0f, 1.0f, 0.0f);
}

void drawCrosshair() {
    glPushMatrix();
    glLoadIdentity();
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    int width = glutGet(GLUT_WINDOW_WIDTH);
    int height = glutGet(GLUT_WINDOW_HEIGHT);
    gluOrtho2D(0, width, 0, height);

    glColor4f(1.0, 0.0, 0.0, 0.3);

    glBegin(GL_LINES);
    glVertex2f(width / 2 - 10, height / 2);
    glVertex2f(width / 2 + 10, height / 2);
    glVertex2f(width / 2, height / 2 - 10);
    glVertex2f(width / 2, height / 2 + 10);
    glEnd();

    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
}

void renderScene() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    updateCamera();
    
    glColor3f(1.0, 1.0, 1.0);
    glPointSize(5.0f);
    glBegin(GL_POINTS);

    for (auto& p : points) {
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

        updateCameraPosition();
    } else if (cam.panning) {
        const float fov = 45.0f;
        float aspectRatio = static_cast<float>(glutGet(GLUT_WINDOW_WIDTH)) / glutGet(GLUT_WINDOW_HEIGHT);
        float baseSensitivity = 0.002f;
        float radianFov = fov * M_PI / 180.0f;
        float effectiveSensitivity = baseSensitivity * (cam.distance * tanf(radianFov / 2) * aspectRatio);

        float dx = static_cast<float>(x - cam.lastX) * -effectiveSensitivity;
        float dy = static_cast<float>(y - cam.lastY) * effectiveSensitivity;

        float rightX = cosf(cam.angleX);
        float rightZ = -sinf(cam.angleX);
        float upY = cosf(cam.angleY);
        float upX = -sinf(cam.angleX) * sinf(cam.angleY);
        float upZ = -cosf(cam.angleX) * sinf(cam.angleY);

        cam.focusX += dx * rightX;
        cam.focusZ += dx * rightZ;
        cam.focusY += dy * upY;
        cam.focusX += dy * upX;
        cam.focusZ += dy * upZ;

        updateCameraPosition();
    }

    cam.lastX = x;
    cam.lastY = y;
    glutPostRedisplay();
}

void mouseButton(int button, int state, int x, int y) {
    if (button == 3 || button == 4) {
        if (state == GLUT_DOWN) {
            cam.distance *= (button == 3) ? 0.9f : 1.1f;
            updateCameraPosition();
            glutPostRedisplay();
        }
    } else {
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

void keyPress(unsigned char key, int x, int y) {
    float stepSize = -0.5f;  // Adjust this value for faster or slower movement
    // Calculate the camera's forward direction vector
    float forwardX = sinf(cam.angleX) * cosf(cam.angleY);
    float forwardY = sinf(cam.angleY);
    float forwardZ = cosf(cam.angleX) * cosf(cam.angleY);

    // Calculate the camera's right direction vector (orthogonal to the forward vector)
    float rightX = cosf(cam.angleX);
    float rightZ = -sinf(cam.angleX);

    switch (key) {
        case 'w':  // Move focus point forward
            cam.focusX += forwardX * stepSize;
            cam.focusY += forwardY * stepSize;
            cam.focusZ += forwardZ * stepSize;
            break;
        case 's':  // Move focus point backward
            cam.focusX -= forwardX * stepSize;
            cam.focusY -= forwardY * stepSize;
            cam.focusZ -= forwardZ * stepSize;
            break;
        case 'a':  // Move focus point left
            cam.focusX += rightX * stepSize;
            cam.focusZ += rightZ * stepSize;
            break;
        case 'd':  // Move focus point right
            cam.focusX -= rightX * stepSize;
            cam.focusZ -= rightZ * stepSize;
            break;
        case ' ':  // Move focus point up (space bar)
            cam.focusY -= stepSize;
            break;
        case 16:  // Move focus point down (Shift key - ASCII code 16)
            cam.focusY += stepSize;
            break;
    }
    updateCameraPosition();  // Update the camera's position
    glutPostRedisplay();     // Redraw the scene with the new camera focus
}


int main(int argc, char **argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowSize(800, 600);
    glutCreateWindow("3D Scatter Plot with Interactive Camera");

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
    glutKeyboardFunc(keyPress);  // Register keyboard handler

    generatePoints(100000);

    glutMainLoop();

    return 0;
}
