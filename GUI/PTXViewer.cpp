/*
Copyright (C) 2010 David Doria, daviddoria@gmail.com

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

// Instantiate and display the GUI

#include <QApplication>

#include "PTXViewerWidget.h"

int main(int argc, char** argv)
{
  QApplication app(argc, argv);

  PTXViewerWidget* ptxViewerWidget;
  if(argc == 2)
    {
    std::cout << "Using constructor with: " << argv[1] << std::endl;
    ptxViewerWidget = new PTXViewerWidget(argv[1]);
    }
  else
    {
    ptxViewerWidget = new PTXViewerWidget;
    }

  ptxViewerWidget->show();

  return app.exec();
}
