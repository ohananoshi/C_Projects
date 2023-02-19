/* Header: cmd_utilities

   A set of possible useful functions for cmd.

   Author: Guilherme Arruda

   GitHub: https://github.com/ohananoshi

   Included Libraries:
    - stdbool.h
    - windows.h
    - wincon.h
*/

#pragma once

#include <stdbool.h>
#include <windows.h>
#include <wincon.h>

//Types: border_shape
//Set of characters to construct boxes
typedef struct{
    short int left_high_corner_chr;
    short int right_high_corner_chr;
    short int left_bottom_corner_chr;
    short int right_botton_corner_chr;
    short int horizontal_character;
    short int vertical_character;
}border_shape;

//Enum: colors
//Set of standard colors for cmd.
enum colors{
    BLACK = 0,
    DARK_BLUE,
    DARK_GREEN,
    DARK_BLUE_GREEN,
    DARK_RED,
    DARK_PURPLE,
    DARK_BROWN,
    DARK_GRAY,
    LIGHT_GRAY,
    LIGHT_BLUE,
    LIGHT_GREEN,
    LIGHT_BLUE_GREEN,
    LIGHT_RED,
    LIGHT_PURPLE,
    LIGHT_YELLOW,
    WHITE
};

border_shape DOUBLE_LINE = {201, 187, 200, 188, 205, 186};

//Function: set_to_fullscreen
//Display the console in fullscreen mode.
bool set_to_fullscreen(){

    //Obter handle da janela atual
    HWND console_window_handle = GetConsoleWindow();
    if(console_window_handle == NULL) return 0;

    //Alterar o estilo da janela para tela cheia
    SetWindowLong(console_window_handle, GWL_STYLE, WS_POPUP | WS_VISIBLE);

    //Redimensionar a janela para as dimensões da tela
    int width = GetSystemMetrics(SM_CXSCREEN);
    int height = GetSystemMetrics(SM_CYSCREEN);
    if(SetWindowPos(console_window_handle, HWND_TOP, 0, 0, width, height, SWP_SHOWWINDOW)) return 0;

    return 1;
}

/*Function: set_text_color
  Change the color of the text (see color_scheme.png)

  To set a background color to selected foreground color,

  just add the background color multiplied by 16 (see example1.png)

  Parameters:
    color - decimal value of color (8bit).
*/
void set_text_color(unsigned short int color){
    SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), color);
}

//Function: move_cursor
//Moves the console cursor for the specified position in coluns and lines.
void move_cursor(unsigned short int x,unsigned short int y){

	COORD coordinates;
	coordinates.X = x;
	coordinates.Y = y;
	SetConsoleCursorPosition(GetStdHandle(STD_OUTPUT_HANDLE), coordinates);

}

//Function: write_text
//Write a text in defined position.
void write_text(char *text, short int x, short int y){
    move_cursor(x,y);
    printf("%s", text);
}

/*Function: draw_box

  Draw a empty static box with some properties using characters ASCII.

  Parameters:

    ascii_dec_symbols - Pointer with a set of characters necessaries for the construction
                        of box.
    pos               - Position of box relative Left high corner of console

*/

void draw_box(border_shape ascii_dec_symbols,
              short int width,
              short int height,
              short int x,
              short int y){

    move_cursor(x, y);

    printf("%c", ascii_dec_symbols.left_high_corner_chr);

    for(short int i = 0; i < (width - 2); i++){
        printf("%c",ascii_dec_symbols.horizontal_character);
    }

    printf("%c", ascii_dec_symbols.right_high_corner_chr);

    for(short int j = 1; j < (height-2); j++){
        move_cursor(x, y + j);
        printf("%c", ascii_dec_symbols.vertical_character);
        move_cursor(x + (width - 1), y + j);
        printf("%c", ascii_dec_symbols.vertical_character);
    }

    move_cursor(x, y + height - 2);
    printf("%c", ascii_dec_symbols.left_bottom_corner_chr);

    for(short int i = 0; i < (width - 2); i++){printf("%c",ascii_dec_symbols.horizontal_character);}

    printf("%c", ascii_dec_symbols.right_botton_corner_chr);
}

/* Function: fill_rectangle

   Print a fill rectangle with size of parameters.

   Parameters:

     ascii_dec_character - decimal value os asacii character
     width               - width of rectangle
     height              - height of rectangle
*/

void fill_rectangle(short int ascii_dec_character, short int width, short int height){
    for (short int i = 0; i < height; i++){
        for(short int j = 0; j < width; j++){
            printf("%c", ascii_dec_character);
        }
    }
}


void print_list(char **chr_list, unsigned short int number_of_itens, unsigned short int x, unsigned short int y,
                unsigned short int space, bool direction){

    short int offset = x;

    if(direction){
        for(short int i = 0; i < number_of_itens; i++){
        write_text(chr_list[i], x, y+i*space);
        }
    }else{
        for(short int i = 0; i < number_of_itens; i++){
        write_text(chr_list[i], x+offset + i*space, y);
        offset = offset + strlen(chr_list[i]);
        }
    }
}
