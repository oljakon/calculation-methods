function lab1()
C = [10, 4, 9, 8, 5;
    9, 3, 5, 7, 8;
    2, 5, 8, 10, 5;
    4, 5, 7, 9, 3;
    8, 7, 10, 9, 6];

debug = true;
%debug = false;
task_max = true;
%task_max = false;

[rows, cols] = size(C);

C1 = C;

if task_max
    % Поиск максимального элемента, вычитание его из элементов
    % матрицы и домножение на -1
    max_cost = 0;
    for i = 1:rows
        for j = 1:cols
            if C(i,j) > max_cost
                max_cost = C(i,j);
            end
        end
    end
    temp_c = C;
    for i = 1:rows
        for j = 1:cols
            C1(i,j) = temp_c(i,j)*(-1) + max_cost;
        end
    end
    
    if debug
        disp("Матрица после первого шага максимизации");
        disp(C1);
    end
end

% Поиск наименьшего элемента столбца и его вычитание из элементов столбца
min_c = min(C1, [], 1);
C1 = C1 - min_c;
if debug
    disp("Матрица после вычитания наименьшего элемента по столбцам");
    disp(C1);
end

% Поиск наименьшего элемента строки и его вычитание из элементов строки
min_r = min(C1, [], 2);
C1 = C1 - min_r;
if debug
    disp("Матрица после вычитания наименьшего элемента по строкам");
    disp(C1);
end

% Первичный поиск 0*
stars = find_stars(C1);

N = size(C1, 1);
quotes = zeros(N);

if debug
    disp("Выделенные 0* для определения CHH");
    debug_matrix(C1, stars, quotes, [], zeros(N));
end
% |CHH|
CNN = nnz(stars);
if debug
    disp("CHH = ");
    disp(CNN);
end
if CNN < cols
    if debug
        disp("СНН нужно улучшать");
    end
    
    cur_row = 0;
    cur_col = 0;
    
    iteration = 1;
    while CNN < rows
        if debug
            disp("Номер итерации: ");
            disp(iteration);
        end

        markedCol = find_marked_cols(stars);
        markedRows = zeros(size(C1,2));
        while true
            % Поиск неотмеченных нулей
            [succ, cur_col, cur_row] = find_unm_zeros(C1, markedCol, markedRows);
            % Если 0 найден
            if succ
                quotes(cur_row,cur_col) = 1;
                if debug
                    disp("Матрица после нахождения неотмеченного нуля");
                    debug_matrix(C1, stars, quotes, markedCol, markedRows);
                end
                % Отметить и продолжить поиск
                [res, col] = check_marked_zero_in_row(cur_row, stars);
                if res == true
                    markedRows(cur_row) = 1;
                    markedCol(col) = 0;
                    if debug
                        disp("Матрица после переопределения выделений строк и столбцов");
                        debug_matrix(C1, stars, quotes, markedCol, markedRows);
                    end
                % Если ноль не найден, построить L-цепочку
                else
                    stars = create_L_chain(stars, quotes, cur_row, cur_col);
                    CNN=CNN+1;
                    markedCol = find_marked_cols(stars);
                    markedRows = zeros(size(C1,2));
                    if debug
                        disp("Матрица 0* после построения L-цепочки: ");
                        disp(stars);
                    end
                    break;
                end
            % Если ноль не найден, вычесть h
            else
                C1 = calc_h(C1, markedRows, markedCol);
                if debug
                    disp("Матрица после вычитания h из столбцов и прибавления к строкам");
                    debug_matrix(C1, stars, quotes, markedCol, markedRows);
                end
            end
        end
    iteration = iteration + 1;
    end
end

% X_opt и f_opt
X_opt = zeros(cols);
for i = 1:rows
    for j = 1:cols
        if stars(i,j)
            X_opt(i,j) = 1;
        end
    end
end
disp("X_opt = ");
disp(X_opt);
f_opt = 0;
for i = 1:rows
    for j = 1:cols
        if stars(i,j)
            f_opt = f_opt + C(i,j)*stars(i,j);
        end
    end
end
disp("f_opt = ");
disp(f_opt);

% Функция первичного поиска 0*
function res = find_stars(C1)
    N = size(C1,1);
    stars = zeros(N);
    for i = 1:N
        for j = 1:N
            if C1(i,j) == 0
                if max(stars(:,j))==0 && max(stars(i,:))==0
                    stars(i,j) = 1;
                    break;
                end
            end
        end
    end
    res = stars;
end

% Функция поиска отмеченных столбцов
function marked_cols = find_marked_cols(stars)
    marked_cols = sum(stars);
end

% Функция поиска неотмеченного нуля в матрице
function [res,col,row] = find_unm_zeros(C1, markedCol, markedRows)
    N = size(C1,1);
    col = 0;
    row = 0;
    res = false;
    for j = 1:N
        for i = 1:N
            if C1(i,j) == 0 && markedCol(j) == 0 && markedRows(i) == 0
                res = true;
                col = j;
                row = i;
                break;
            end
        end
    end
end

% Функция, которая определяет, есть ли в строке матрицы 0*
function [res, col] = check_marked_zero_in_row(row, stars)
    res = false;
    col = 0;
    for j = 1:size(stars, 2)
        if stars(row, j) == 1
            res = true;
            col = j;
            break
        end
    end
end

% Функция поиска 0* в столбце
function [res, row] = check_marked_zero_in_column(col,stars)
    [res,row] = check_marked_zero_in_row(col,stars');
end

% Функция пересчета матрицы на основе поиска h и его вычитания из невыделенных столбцов
function matrix = calc_h(C1, markedRows, markedCol)
    matrix = C1;
    N = size(C1,1);
    min_el = -1;
    
    % Поиск минимального элемента из невыделенных
    for j = 1:N
        for i = 1:N
            if markedCol(j) == 0 && markedRows(i) == 0
                if matrix(i, j) < min_el || min_el == -1
                    min_el = matrix(i, j);
                end
            end
        end
    end
    
    % Вычитание минимального элемента из невыделенных
    for j = 1:N
        for i = 1:N
            if markedCol(j) == 0 && markedRows(i) == 0
                matrix(i,j) = matrix(i,j) - min_el;
            end
            if markedCol(j) == 1 && markedRows(i) == 1
                matrix(i,j) = matrix(i,j) + min_el;
            end
        end
    end
end

% Функция построения L-цепочки
function created_L = create_L_chain(stars, quotes, row, col)
    res = true;
    cur_col = col;
    cur_row = row;
    created_L = stars;
    while res
        created_L(cur_row,cur_col) = 1;
        [res,cur_row] = check_marked_zero_in_column(cur_col,stars);
        if(res == true)
            stars(cur_row, cur_col) = 0;
            created_L(cur_row, cur_col) = 0;
            [res, cur_col] = check_marked_zero_in_row(cur_row,quotes);
        end
    end
end

% Функция вывода матрицы с 0*, 0'
function debug_matrix(C1, stars, quotes, markedCol, markedRows)
    for i = 1:size(C1,1)
        for j = 1:size(C1,2)
            if stars(i,j) == 1
                fprintf("%d* \t", C1(i,j));
            else
                if quotes(i,j) == 1
                    fprintf("%d' \t", C1(i,j));
                else
                    fprintf("%d \t", C1(i,j));
                end
            end
        end
        if markedRows(i)
            fprintf(' +\n');
        else
            fprintf(' \n');
        end 
    end
    for i = 1:length(markedCol)
        if markedCol(i)
            fprintf("+ \t");
        else
            fprintf(" \t");
        end
    end
    fprintf("\n");
end

end