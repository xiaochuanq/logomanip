#include <iostream>
#include <getopt.h>
#include "cv.h"
#include "highgui.h"

using namespace std;
using namespace cv;

int main( int argc, char** argv )
{
  if( argc < 3 ){
    cout << "Usage: " << argv[0] <<" image logo clipart" << endl;
    cout << "image:   the image to be tested" << endl;
    cout << "logo:    the logo image" << endl;
    cout << "clipart: the clip art image" <<endl;
  }

  Mat testimg = imread( argv[1] );
  Mat logoimg = imread( argv[2] );
  Mat clipimg = imread( argv[3] );

  Mat testimg_gray, logoimg_gray;
  cvtColor( testimage, testimg_gray, CV_RGB2GRAY );
  cvtColor( logoimage, logoimg_gray, CV_RGB2GRAY );
}
