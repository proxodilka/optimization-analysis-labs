function pretty_print(x)
    printf("[")
    for i=1:length(x)
        printf("%e", x(i))
        if i ~= length(x) then printf(", ") end
    end
    printf("]")
endfunction

function [res]=to_list(arr)
    res = list()
    for i=1:length(arr)
        res($ + 1) = arr(i)
    end
endfunction

function [res]=simple_proj(point, square)
    [x, y] = to_list(point)(:)
    [x_diff, y_diff] = list(0, 0)(:)

    if x < square(1) then
        x_diff = square(1) - x
    elseif x > square(2) then
        x_diff = square(2) - x
    end
    if y < square(3) then
        y_diff = square(3) - y
    elseif y > square(4) then
        y_diff = square(4) - y
    end

    if size(point)(1) == 1 && size(point)(2) != 1
        point = point'
    end

    res = point + [x_diff, y_diff]'
endfunction

function [res]=f1(point)
    [x, y] = to_list(point)(:)
    res = 100 * (x^2 - y)^2 + (1-x)^2
endfunction

function [res]=f2(point)
    [x, y] = to_list(point)(:)
    res = 100 * (y - x^3)^2 + (1-x)^2
endfunction

function [res]=f2_grad(point)
    [x, y] = to_list(point)(:)
    res = [
        2 * (300*x^5 - 300*x^2*y + x - 1),
        200 * (y - x^3)
    ]'
endfunction

function [res]=f2_proj(point)
    res = simple_proj(point, list(-1.2, 1, -1, 1))'
endfunction

function [res]=f4(point)
    [x1, x2, x3, x4] = to_list(point)(:)
    res = (x1 + 10*x2)^2 + 5*(x3-x4) ^ 2 + (x2 - 2*x3)^4 + 10*(x1 - x4)^4
endfunction

function [res]=search(f, x, i, h, eps, div_to)
    base = f(x)
    dir_vector = zeros(x)
    dir_vector(i) = 1
    base_h = h
    
    test_point = x + dir_vector * h
    while (%T)
        if (h <= eps)
            break
        end
        h = h / div_to
        test_point = x + dir_vector * h
        if f(test_point) < base
            break
        end
        test_point = x - dir_vector * h
        if f(test_point) < base
            break
        end
    end
    
    if f(test_point) < base then
        res = test_point
    else
        res = x
    end
endfunction

function [res]=hooke_jeeves(f, x0, max_step_scale, scale_reducer, stop_condition, projection, verbose_level, max_iter)
    // Минимизирует целевую функцию методом Хука — Дживса.
    //
    // Аргументы
    // ---------
    // f : callable(vector) -> float
    //    Целевая функция для минимизации.
    // x0 : vector
    //    Начальное приближение.
    // max_step_scale : float, по умолчанию: 1
    //    Максимально возможный коэффициент шага метода.
    // scale_reducer : float, по умолчанию: 2
    //    Коэффициент уменьшения шага метода в случае "перепрыгивания" точки минимума.
    // stop_condition : float, callable(f, vector, vector) -> bool, по умолчанию: 10e-7
    //    Если передано число, параметр интерпретируется как точность метода, критерий останова при
    //    этом будет следующим:
    //      `` abs(f(xk) - f(xk_1)) < stop_condition ``
    //    Если передана функция, параметр интерпретируется как критерий останова. На каждой итерации
    //    функция получает на вход 3 параметра:
    //      1. Целевую функцию f.
    //      2. Значение x на текущем шаге.
    //      3. Значение х на предыдущем шаге.
    //    Если функция возвращает True, метод прекращает свою работу.
    // projection : callable(vector) -> vector, опционально
    //    Проекция на допустимое множество решения. Метод гарантирует, что все рассматриваемые
    //    точки будут находится внутри этого множества. Можно использовать для решения задач 
    //    условного экстремума.
    // verbose_level : int, по умолчанию: 1
    //    Регулирует количество информации информации о работе метода, выводимой в консоль:
    //          0: Не выводить информацию в консоль.
    //        >=1: Вывести в консоль результат работы метода.
    //        >=2: Выводить в консоль информацию о текущем решении на каждом шаге.
    // max_iter : int, по умолчанию: 1000
    //    Максимальное число итераций метода.
    //
    // Возвращает
    // ----------
    // vector
    //    Найденная точка минимума функции.
    //
    // Примечание
    // ----------
    // Ссылка на теоретическое описание работы метода: https://intuit.ru/studies/courses/1020/188/lecture/4931

    // ========================= Обработка значений по умолчанию =========================
    if ~exists("max_step_scale", "local") then max_step_scale = 1 end
    if ~exists("scale_reducer", "local") then scale_reducer = 2 end
    if ~exists("verbose_level", "local") then verbose_level = 1 end
    if ~exists("max_iter", "local") then max_iter = 1000 end
    if ~exists("eps", "local") then eps = 10e-7 end
    if ~exists("stop_condition", "local") then
        function [_res] = stop_condition(f, xk, xk_1)
            _res = abs(norm(xk - xk_1)) < eps
        endfunction
    end
    if ~exists("projection", "local") then
        function [_res] = projection(point)
            _res = point
        endfunction
    end

    // ========================= Начало метода =========================
    better_x = x0
    for j=1:max_iter
        scale = max_step_scale
        better_x = x0
        prev_x = x0
        for i=1:length(x0)
            better_x = projection(search(f, better_x, i, scale, eps, scale_reducer))
        end
        direction = better_x - x0
        i = 1
        while (scale >= eps)
            while f(projection(x0 + scale * direction * i)) < f(projection(x0 + scale * direction * (i - 1)))
                i = i + 1
            end
            x0 = projection(x0 + scale * direction * (i - 1))
            i = 1
            scale = scale / scale_reducer
        end
        x0 = projection(x0 + scale * direction * (i - 1))

        if verbose_level > 1 then
            printf("Step %i: xk = ", j)
            pretty_print(x0)
            printf("| F(xk) = %e\n", f(x0))
        end

        if stop_condition(f, x0, prev_x)
            break
        end
    end

    if verbose_level > 0 then
        printf("Found solution in %i iterations:\n", j)
        printf("\tResult point: ")
        pretty_print(x0)
        printf("\n\tFunction value: %e\n", f(x0))
    end
    res = x0
endfunction

function [res]=golden_section_search(fn, sect, eps, max_iter)
    // Минимизирует целевую функцию одного переменного методом золотого сечения.
    //
    // Аргументы
    // ---------
    // f : callable(float) -> float
    //    Целевая функция для минимизации.
    // sect : list
    //    Отрезок для поиска минимизирующей точки.
    // eps : float
    //    Погрешность метода.
    // max_iter : int
    //    Максимальное число итераций метода.
    //
    // Возвращает
    // ----------
    // float
    //    Найденная точка минимума функции.
    //
    // Примечание
    // ----------
    // Ссылка на теоретическое описание работы метода: http://www.machinelearning.ru/wiki/index.php?title=%D0%9C%D0%B5%D1%82%D0%BE%D0%B4_%D0%B7%D0%BE%D0%BB%D0%BE%D1%82%D0%BE%D0%B3%D0%BE_%D1%81%D0%B5%D1%87%D0%B5%D0%BD%D0%B8%D1%8F._%D0%A1%D0%B8%D0%BC%D0%BC%D0%B5%D1%82%D1%80%D0%B8%D1%87%D0%BD%D1%8B%D0%B5_%D0%BC%D0%B5%D1%82%D0%BE%D0%B4%D1%8B1
    [a, b] = to_list(sect)(:)
    
    ff = (1 + sqrt(5)) / 2
    x1 = b - (b - a) / ff
    x2 = a + (b - a) / ff
    for i=1:max_iter
        if fn(x1) > fn(x2) then
            a = x1
            x1 = x2
            x2 = b - (x1 - a)
        end
        if fn(x1) < fn(x2) then
            b = x2
            x2 = x1
            x1 = a + (b - x2)
        end
        if (fn(x1) == fn(x2)) then
            a = x1
            b = x1
        end
        
        res = (a + b) / 2
        
        err = (b - a) / 2
        if err < eps then break end
    end
endfunction

function [x0]=grad_descent(f, x0, max_step, grad_fn, projection, stop_condition, verbose_level, max_iter)
    // Минимизирует целевую функцию градиентным методом наискорейшего спуска.
    //
    // Аргументы
    // ---------
    // f : callable(vector) -> float
    //    Целевая функция для минимизации.
    // x0 : vector
    //    Начальное приближение.
    // max_step : float, по умолчанию: 100
    //    Максимально возможный коэффициент шага метода.
    // grad_fn : callable(vector) -> vector, опционально
    //    Градиент целевой функции, если не передано, будет вычислено автоматически.
    // projection : callable(vector) -> vector, опционально
    //    Проекция на допустимое множество решения. Метод гарантирует, что все рассматриваемые
    //    точки будут находится внутри этого множества. Можно использовать для решения задач 
    //    условного экстремума.
    // stop_condition : float, callable(f, vector, vector) -> bool, по умолчанию: 10e-7
    //    Если передано число, параметр интерпретируется как точность метода, критерий останова при
    //    этом будет следующим:
    //      `` abs(f(xk) - f(xk_1)) < stop_condition ``
    //    Если передана функция, параметр интерпретируется как критерий останова. На каждой итерации
    //    функция получает на вход 3 параметра:
    //      1. Целевую функцию f.
    //      2. Значение x на текущем шаге.
    //      3. Значение х на предыдущем шаге.
    //    Если функция возвращает True, метод прекращает свою работу.
    // verbose_level : int, по умолчанию: 1
    //    Регулирует количество информации информации о работе метода, выводимой в консоль:
    //          0: Не выводить информацию в консоль.
    //        >=1: Вывести в консоль результат работы метода.
    //        >=2: Выводить в консоль информацию о текущем решении на каждом шаге.
    // max_iter : int, по умолчанию: 1000
    //    Максимальное число итераций метода.
    //
    // Возвращает
    // ----------
    // vector
    //    Найденная точка минимума функции.
    //
    // Примечание
    // ----------
    // Ссылка на теоретическое описание работы метода: http://www.machinelearning.ru/wiki/index.php?title=%D0%9C%D0%B5%D1%82%D0%BE%D0%B4_%D0%B3%D1%80%D0%B0%D0%B4%D0%B8%D0%B5%D0%BD%D1%82%D0%BD%D0%BE%D0%B3%D0%BE_%D1%81%D0%BF%D1%83%D1%81%D0%BA%D0%B0#.D0.9C.D0.B5.D1.82.D0.BE.D0.B4_.D0.BD.D0.B0.D0.B8.D1.81.D0.BA.D0.BE.D1.80.D0.B5.D0.B9.D1.88.D0.B5.D0.B3.D0.BE_.D1.81.D0.BF.D1.83.D1.81.D0.BA.D0.B0_2
    
    // ========================= Обработка значений по умолчанию =========================
    if ~exists("max_step", "local") then max_step = 100 end
    if ~exists("eps", "local") then eps=10e-7 end
    if ~exists("verbose_level") then verbose_level = 1 end
    if ~exists("max_iter", "local") then max_iter = 1000 end
    if ~exists("grad_fn", "local") then 
        function [_res] = grad_fn(x)
            _res = numderivative(f, x)
        endfunction
    end
    if ~exists("stop_condition", "local") then
        function [_res] = stop_condition(f, xk, xk_1)
            _res = abs(norm(xk - xk_1)) < eps
        endfunction
    end
    if ~exists("projection", "local") then
        function [_res] = projection(x)
            _res = x'
        endfunction
    end

    // ========================= Начало метода =========================
    for i=1:max_iter
        grad = grad_fn(x0)'
        function [_res]=get_diff(p)
            _res = f(x0 - p * grad)
        endfunction
        a = golden_section_search(get_diff, list(0, max_step), 10e-5, 1000)
        x_prev = projection(x0)'
        x0 = projection(x0 - a * grad)'
        
        if verbose_level > 1 then
            printf("Step %i: xk = ", i)
            pretty_print(x0)
            printf("| F(xk) = %e\n", f(x0))
        end
        
        if stop_condition(f, x0, x_prev) then break end
    end
    
    if verbose_level > 0 then
        printf("Found solution in %i iterations:\n", i)
        printf("\tResult point: ")
        pretty_print(x0)
        printf("\n\tFunction value: %e\n", f(x0))
    end

endfunction

function [res]=nelder_mead(f, simplex, a, b, c, projection, stop_condition, verbose_level, max_iter)
    // Минимизирует целевую функцию методом Нелдера-Мида.
    //
    // Аргументы
    // ---------
    // f : callable(vector) -> float
    //    Целевая функция для минимизации.
    // simplex : list of vectors
    //    Начальный симплекс.
    // a : float, по умолчанию: 1
    //    Коэффициент отражения.
    // b : float, по умолчанию: 0.5
    //    Коэффициент сжатия.
    // c : float, по умолчанию: 2
    //    Коэффициент растяжения.
    // projection : callable(vector) -> vector, опционально
    //    Проекция на допустимое множество решения. Метод гарантирует, что все рассматриваемые
    //    точки будут находится внутри этого множества. Можно использовать для решения задач 
    //    условного экстремума.
    // stop_condition : float, callable(f, simplex, simplex_prev) -> bool, по умолчанию: 10e-7
    //    Если передано число, параметр интерпретируется как точность метода, критерий останова при
    //    этом будет следующим:
    //      `` norm(simplex -simplex_prev) < stop_condition ``
    //    Если передана функция, параметр интерпретируется как критерий останова. На каждой итерации
    //    функция получает на вход 2 параметра:
    //      1. Целевую функцию f.
    //      2. Симплекс на текущем шаге.
    //      3. Симплекс на предыдущем шаге.
    //    Если функция возвращает True, метод прекращает свою работу.
    // verbose_level : int, по умолчанию: 1
    //    Регулирует количество информации информации о работе метода, выводимой в консоль:
    //          0: Не выводить информацию в консоль.
    //        >=1: Вывести в консоль результат работы метода.
    //        >=2: Выводить в консоль информацию о текущем решении на каждом шаге.
    //        >=3: Выводить в консоль информацию о текущем симплексе на каждом шаге.
    // max_iter : int, по умолчанию: 1000
    //    Максимальное число итераций метода.
    //
    // Возвращает
    // ----------
    // vector
    //    Найденная точка минимума функции.
    //
    // Примечание
    // ----------
    // Ссылка на теоретическое описание работы метода: http://www.machinelearning.ru/wiki/index.php?title=%D0%9C%D0%B5%D1%82%D0%BE%D0%B4_%D0%9D%D0%B5%D0%BB%D0%B4%D0%B5%D1%80%D0%B0-%D0%9C%D0%B8%D0%B4%D0%B0

    // ========================= Обработка значений по умолчанию =========================
    if ~exists("a", "local") then a = 1 end
    if ~exists("b", "local") then b = 0.5 end
    if ~exists("c", "local") then c = 2 end
    if ~exists("verbose_level", "local") then verbose_level = 1 end
    if ~exists("max_iter", "local") then max_iter = 1000 end
    if ~exists("eps", "local") then eps = 10e-7 end
    if ~exists("stop_condition", "local") then
        function [_res] = stop_condition(f, simplex_k, simplex_k_1)
            max_dist = -%inf
            indices = [1, length(simplex_k) - 1, length(simplex_k)]
            for j=1:length(indices)
                i = indices(j)
                cur_dist = abs(norm(simplex_k(i) - simplex_k_1(i)))
                if cur_dist > max_dist then max_dist = cur_dist end
            end
            _res = abs(max_dist) < eps
        endfunction
    end
    if ~exists("projection", "local") then
        function [_res] = projection(point)
            _res = point
        endfunction
    end

    // ========================= Вспомогательные функции =========================
    function [_res] = sort_simplex(simplex)
        simplex_values = zeros(length(simplex))
        for i=1:length(simplex)
            simplex_values(i) = f(simplex(i))
        end
        [_, indices] = gsort(simplex_values, "g", "i")
        _res = list(simplex(indices))
    endfunction

    function print_point(point, i, prepend)
        if ~exists("prepend", "local") then prepend = "" end
        if type(i) ~= 10 then i = string(i) end

        printf(prepend + i + ": ")
        pretty_print(point)
        printf(" | %f\n", f(point))
    endfunction

    function print_simplex(simplex, verbose, prepend)
        if ~exists("verbose", "local") then verbose = %F end
        if ~exists("prepend", "local") then prepend = "" end

        if ~verbose
            for i=1:length(simplex)
                print_point(simplex(i), i, prepend)
            end
        else
            print_point(simplex(1), "L", prepend)
            if length(simplex) > 3 then printf(prepend + "...\n") end
            print_point(simplex($-1), "G", prepend)
            print_point(simplex($), "H", prepend)
        end
    endfunction

    // ========================= Начало метода =========================
    for j=1:max_iter
        prev_simplex = simplex
        simplex = sort_simplex(simplex)
        accum = zeros(simplex(1))

        for i=1:length(simplex) - 1
            accum = accum + simplex(i)
        end

        [xh, xg, xl] = list(simplex($), simplex($ - 1), simplex(1))(:)
        x_prev = xl

        xc = accum ./ (length(simplex) - 1)
        xr = projection((1 + a) * xc - a * xh)

        should_squeeze = %F
        if f(xr) < f(xl) then
            xs = projection((1 - c) * xc + c * xr)
            if f(xs) < f(xl) then
                xh = xs
            else
                xh = xr
            end
        elseif f(xl) < f(xr) & f(xr) < f(xg) then
            xh = xr
        elseif f(xh) > f(xr) & f(xr) > f(xg) then
            [xr, xh] = list(xh, xr)(:)
            should_squeeze = %T
        else
            should_squeeze = %T
        end

        if should_squeeze
            xs = projection(b * xh + (1 - b) * xc)
            if f(xs) < f(xh)
                xh = xs
                simplex($) = xh
            else
                for i=1:length(simplex)
                    simplex(i) = projection(xl + (simplex(i) - xl) ./ 2)
                end
            end
        end

        simplex($) = xh
        
        if f(simplex(1)) < f(simplex($))
            res = simplex(1)
        else
            res = simplex($)
        end

        if verbose_level > 1 then
            printf("Step %i: xk = ", j)
            pretty_print(res)
            printf("| F(xk) = %e\n", f(res))
            if verbose_level > 2 then
                printf("Current simplex:\n")
                print_simplex(sort_simplex(simplex), verbose=%T, prepend="\t")
                printf("\n")
            end
        end

        if stop_condition(f, simplex, prev_simplex) then break end
    end

    if verbose_level > 0 then
        printf("Found solution in %i iterations:\n", j)
        printf("\tResult point: ")
        pretty_print(res)
        printf("\n\tFunction value: %e\n", f(res))
    end
endfunction

eps = 10e-7
max_iter = 1000

function [res]=solid_stop_condition(f, xk, xk_1)
    res = abs(f(xk) - f(xk_1)) < eps && abs(norm(xk - xk_1)) < eps
endfunction

printf("\n============================ Метод Хука-Дживса ============================\n")

printf("\n--------------- Problem 1 ---------------\n")
printf("100(x^2 - y)^2 + (1-x)^2 -> min | X0 = (-1.2, 1) | eps = 10e-7\n")
hooke_jeeves(f1, [-1.2, 1], max_step_scale=0.5, scale_reducer=2, stop_condition=solid_stop_condition)

printf("\n--------------- Problem 2 ---------------\n")
printf("100(y - x^3)^2 + (1-x)^2 -> min | X0 = (-1.2, -1) | eps = 10e-7\n")
hooke_jeeves(f2, [-1.2, -1], scale_reducer=1.2, stop_condition=solid_stop_condition)

printf("\n--------------- Problem 3 ---------------\n")
printf("100(y - x^3)^2 + (1-x)^2 -> min | X0 = (-1.2, -1) | eps = 10e-7 | x ∈ [-1.2, 1] | y ∈ [-1, 1]\n")
hooke_jeeves(f2, [-1.2, -1], scale_reducer=1.2, stop_condition=solid_stop_condition, projection=f2_proj)

printf("\n--------------- Problem 4 ---------------\n")
printf("(x1 + 10x2)^2 + 5(x3 - x4)^2 + (x2 - 2x3)^4 + 10(x1 - x4)^4 | X0 = (3, -1, 0, 1) | eps = 10e-7\n")
hooke_jeeves(f4, [3, -1, 0, 1], max_step_scale=0.5, scale_reducer=1.5, stop_condition=solid_stop_condition)


printf("\n============================ Метод Нелдера-Мида ============================\n")

printf("\n--------------- Problem 1 ---------------\n")
printf("100(x^2 - y)^2 + (1-x)^2 -> min | X0 = [(-1.2, 1), (3, -2), (0, 0)] | eps = 10e-7\n")
nelder_mead(f1, simplex=list([-1.2, 1], [3, -2], [0, 0]))

printf("\n--------------- Problem 2 ---------------\n")
printf("100(y - x^3)^2 + (1-x)^2 -> min | X0 = [(-1.2, -1), (-3, 2), (0, 0)] | eps = 10e-7\n")
nelder_mead(f2, simplex=list([-1.2, -1], [-3, 2], [0, 0]))

printf("\n--------------- Problem 3 ---------------\n")
printf("100(y - x^3)^2 + (1-x)^2 -> min | X0 = [(-1.2, -1), (-1.2, 1), (0, 0)] | eps = 10e-7 | x ∈ [-1.2, 1] | y ∈ [-1, 1]\n")
nelder_mead(f2, simplex=list([-1.2, -1], [-1.2, 1], [0, 0]), projection=f2_proj)

printf("\n--------------- Problem 4 ---------------\n")
printf("(x1 + 10x2)^2 + 5(x3 - x4)^2 + (x2 - 2x3)^4 + 10(x1 - x4)^4 | X0 = [(3, -1, 0, 1), (-3, 1, 0, -1), (0, 0, 1, 2), (1, 2, 0, 0), (1, 2, 3, 4)] | eps = 10e-7\n")
nelder_mead(f4, simplex=list([3, -1, 0, 1], [-3, 1, 0, -1], [0, 0, 1, 2], [1, 2, 0, 0], [1, 2, 3, 4]))


printf("\n============================ Градиентный метод наискорейшего спуска ============================\n")

printf("\n--------------- Problem 1 ---------------\n")
printf("100(x^2 - y)^2 + (1-x)^2 -> min | X0 = (-1.2, 1) | eps = 10e-7\n")
grad_descent(f1, x0=[-1.2, 1]', max_step=200, stop_condition=solid_stop_condition)

printf("\n--------------- Problem 2 ---------------\n")
printf("100(y - x^3)^2 + (1-x)^2 -> min | X0 = (-1.2, -1) | eps = 10e-7\n")
grad_descent(f2, x0=[-1.2, -1]', stop_condition=solid_stop_condition)

printf("\n--------------- Problem 3 ---------------\n")
printf("100(y - x^3)^2 + (1-x)^2 -> min | X0 = (-1.2, -1) | eps = 10e-7 | x ∈ [-1.2, 1] | y ∈ [-1, 1]\n")
grad_descent(f2, x0=[-1.2, -1]', max_step=200, grad_fn = f2_grad, projection=f2_proj, stop_condition=solid_stop_condition)

printf("\n--------------- Problem 4 ---------------\n")
printf("(x1 + 10x2)^2 + 5(x3 - x4)^2 + (x2 - 2x3)^4 + 10(x1 - x4)^4 | X0 = (3, -1, 0, 1) | eps = 10e-7\n")
grad_descent(f4, x0=[3, -1, 0, 1]', stop_condition=solid_stop_condition, max_iter=10000)
