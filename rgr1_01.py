import math
import itertools
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

# --- Геометрические вспомогательные функции ---

def cross_product(p1, p2, p3):
    """
    Вычисляет векторное произведение (p2-p1) x (p3-p1).
    Если результат > 0, то p3 находится слева от вектора p1->p2 (против часовой стрелки).
    Если результат < 0, то p3 находится справа от вектора p1->p2 (по часовой стрелке).
    Если результат = 0, то точки коллинеарны.
    """
    return (p2[0] - p1[0]) * (p3[1] - p1[1]) - (p2[1] - p1[1]) * (p3[0] - p1[0])

def distance_sq(p1, p2):
    """Квадрат расстояния между двумя точками."""
    return (p1[0] - p2[0])**2 + (p1[1] - p2[1])**2

def convex_hull_graham_scan(points):
    """
    Строит минимальную выпуклую оболочку для набора точек с помощью сканирования Грэхема.
    Возвращает список вершин оболочки в порядке против часовой стрелки.
    """
    if len(points) < 3:
        # Недостаточно точек для формирования многоугольника.
        # Возвращаем отсортированные точки, чтобы избежать ошибок при попытке доступа к элементам.
        # Для данной задачи, если < 3, это не будет четырехугольник.
        return sorted(points, key=lambda p: (p[0], p[1]))


    # 1. Найти опорную точку p0 (самая нижняя, затем самая левая).
    p0 = min(points, key=lambda p: (p[1], p[0]))

    # 2. Отсортировать остальные точки по полярному углу относительно p0.
    # Если углы равны, берем более дальнюю точку.
    def sort_key(p):
        angle = math.atan2(p[1] - p0[1], p[0] - p0[0])
        # Для коллинеарных точек, более дальняя точка должна идти позже
        return (angle, -distance_sq(p0, p))

    sorted_points = sorted([p for p in points if p != p0], key=sort_key)
    
    # Добавляем p0 в начало отсортированного списка
    # (или используем points отсортированные относительно p0, включая p0)
    # Для простоты, создадим новый список с p0 в начале
    
    # Отсортировать все точки (кроме p0) по полярному углу относительно p0
    # Если две точки имеют одинаковый угол, отбрасываем ту, что ближе к p0
    # (или обрабатываем коллинеарные точки аккуратно в основном цикле сканирования)
    
    # Для простоты, используем стандартный подход:
    # Сначала p0, затем остальные, отсортированные по углу.
    # Если есть точки с одинаковым углом, берем самую дальнюю.
    
    # Пересортируем все точки (включая p0, если она не первая после сортировки по y,x)
    # относительно новой p0.
    
    # Создаем список точек для сортировки (все, кроме p0)
    other_points = [p for p in points if p != p0]
    
    # Сортируем other_points
    # При одинаковых углах, берем самую дальнюю точку.
    # math.atan2 возвращает угол в [-pi, pi].
    # Для обработки коллинеарных точек, если углы одинаковы,
    # сначала должна идти более близкая точка, чтобы дальняя ее "перекрыла" при необходимости.
    # Либо, если углы одинаковы, берем самую дальнюю (стандартный Грэхем).
    # Для простоты, если углы равны, сортируем по расстоянию (ближайшая первая).
    # Сканирование Грэхема само разберется с коллинеарными точками на пути.
    
    def angle_sort_key(point):
        angle = math.atan2(point[1] - p0[1], point[0] - p0[0])
        # Приоритет по углу, затем по расстоянию (для обработки коллинеарных)
        return (angle, distance_sq(p0, point))

    sorted_angular_points = sorted(other_points, key=angle_sort_key)


    hull = [p0, sorted_angular_points[0]] # Инициализируем стек первыми двумя точками

    for i in range(1, len(sorted_angular_points)):
        pt = sorted_angular_points[i]
        # Пока поворот не левый (т.е. правый или коллинеарный)
        while len(hull) >= 2 and cross_product(hull[-2], hull[-1], pt) <= 0:
            hull.pop() # Удаляем последнюю точку из стека
        hull.append(pt)
    
    return hull


def polygon_area(vertices):
    """
    Вычисляет площадь многоугольника по формуле шнурков (Shoelace formula).
    Вершины должны быть упорядочены (по или против часовой стрелки).
    """
    n = len(vertices)
    if n < 3:
        return 0.0  # Не многоугольник
    
    area = 0.0
    for i in range(n):
        j = (i + 1) % n  # Следующая вершина, замыкающаяся на первую
        area += vertices[i][0] * vertices[j][1]
        area -= vertices[j][0] * vertices[i][1]
    return abs(area) / 2.0

def is_inside_convex_polygon(point, polygon_vertices):
    """
    Проверяет, находится ли точка внутри или на границе выпуклого многоугольника.
    Вершины многоугольника должны быть упорядочены (например, против часовой стрелки).
    """
    n = len(polygon_vertices)
    if n < 3:
        return False # Не многоугольник

    # Знак векторного произведения должен быть одинаковым для всех ребер.
    # Для многоугольника против часовой стрелки, точка слева от каждого ребра (или на нем).
    # cross_product(p1, p2, point) > 0 => точка слева от p1->p2
    # cross_product(p1, p2, point) == 0 => точка коллинеарна p1,p2
    
    expected_sign = 0 # 0 - не определен, 1 - положительный, -1 - отрицательный

    for i in range(n):
        p1 = polygon_vertices[i]
        p2 = polygon_vertices[(i + 1) % n]
        
        cp = cross_product(p1, p2, point)
        
        current_sign = 0
        # Используем небольшую эпсилон для сравнения с нулем из-за погрешностей float
        epsilon = 1e-9 
        if cp > epsilon:
            current_sign = 1
        elif cp < -epsilon:
            current_sign = -1
        # else current_sign = 0 (коллинеарна)

        if current_sign == 0: # Точка на прямой, содержащей ребро
            # Для выпуклого многоугольника, если точка на прямой ребра,
            # она считается "внутри" или "на границе".
            # Дополнительно можно проверить, лежит ли она на отрезке, но для выпуклого
            # достаточно проверки со всеми ребрами.
            continue

        if expected_sign == 0:
            expected_sign = current_sign
        elif expected_sign != current_sign:
            return False # Точка с другой стороны от какого-то ребра

    # Если expected_sign так и не был установлен (все cp были ~0),
    # значит точка коллинеарна всем ребрам.
    # Это возможно, если многоугольник вырожден (линия) и точка на ней,
    # или если точка является одной из вершин.
    # В таких случаях считаем, что точка внутри/на границе.
    return True

# --- Основная логика ---

def find_min_area_covering_quadrilateral(set1, set2):
    """
    Находит 4 точки из set1, образующие выпуклый четырехугольник минимальной площади,
    который покрывает все точки из set2.
    """
    min_area = float('inf')
    optimal_quad_vertices = None

    # Рассматриваем все комбинации из 4 различных точек из set1
    if len(set1) < 4:
        return None, float('inf') # Недостаточно точек в set1

    for p_combination in itertools.combinations(set1, 4):
        # p_combination - это кортеж из 4 точек
        
        # 1. Построить выпуклую оболочку для этих 4 точек.
        # Это даст нам упорядоченные вершины кандидата в четырехугольники.
        current_points_list = list(p_combination)
        hull_vertices = convex_hull_graham_scan(current_points_list)

        # 2. Убедиться, что оболочка является четырехугольником (т.е. имеет 4 вершины).
        # Это означает, что все 4 выбранные точки являются вершинами их выпуклой оболочки,
        # и они образуют выпуклый четырехугольник.
        if len(hull_vertices) != 4:
            continue # Эта комбинация не образует выпуклый четырехугольник

        # 3. Вычислить площадь полученного четырехугольника
        current_area = polygon_area(hull_vertices)

        # 4. Проверить, покрывает ли этот четырехугольник все точки из set2
        all_points_in_set2_covered = True
        if not set2: # Если set2 пуст, условие покрытия выполнено
            pass
        else:
            for point_s2 in set2:
                if not is_inside_convex_polygon(point_s2, hull_vertices):
                    all_points_in_set2_covered = False
                    break
        
        # 5. Если покрывает и площадь меньше минимальной, обновить результат
        if all_points_in_set2_covered:
            if current_area < min_area:
                min_area = current_area
                optimal_quad_vertices = hull_vertices
                
    return optimal_quad_vertices, min_area

# --- Данные и выполнение ---

if __name__ == '__main__':
    # Пример данных (можно заменить на свои)
    set1_points = [(0,0), (5,0), (5,5), (0,5), (2,2), (1,3), (3,1)]
    set2_points = [(1,1), (2,3), (3,2)]

    # set1_points = [
    #     (0, 0), (10, 0), (10, 10), (0, 10),  # Большой квадрат
    #     (1, 5), (5, 1), (9, 5), (5, 9),      # Внутренние точки
    #     (2,2), (8,2), (8,8), (2,8)           # Еще точки для выбора
    # ]
    # set2_points = [
    #     (3, 3), (7, 7), (4, 6), (6, 4)
    # ]
    
    # Другой пример
    # set1_points = [(0,0), (1,3), (4,4), (5,1), (2,-1)]
    # set2_points = [(2,2), (3,1)]

    # Пример, где нет решения
    # set1_points = [(0,0), (1,0), (0,1), (1,1)]
    # set2_points = [(5,5)]


    optimal_quad, min_area_val = find_min_area_covering_quadrilateral(set1_points, set2_points)

    # --- Графическое представление ---
    plt.figure(figsize=(8, 8))

    # Отобразить точки из set1
    if set1_points:
        s1_x, s1_y = zip(*set1_points)
        plt.scatter(s1_x, s1_y, c='blue', marker='o', label='Множество 1', s=50, zorder=2)

    # Отобразить точки из set2
    if set2_points:
        s2_x, s2_y = zip(*set2_points)
        plt.scatter(s2_x, s2_y, c='red', marker='x', label='Множество 2', s=70, zorder=3)

    # Отобразить оптимальный четырехугольник
    if optimal_quad:
        print(f"Найден оптимальный четырехугольник с вершинами: {optimal_quad}")
        print(f"Минимальная площадь: {min_area_val:.2f}")
        polygon_patch = Polygon(optimal_quad, closed=True, edgecolor='green', facecolor='lightgreen', alpha=0.5, linewidth=2, zorder=1)
        plt.gca().add_patch(polygon_patch)
    else:
        print("Не удалось найти четырехугольник, удовлетворяющий условиям.")

    plt.xlabel("X координата")
    plt.ylabel("Y координата")
    plt.title("Минимальный покрывающий четырехугольник")
    plt.legend()
    plt.grid(True)
    plt.axis('equal') # Для корректного отображения геометрии
    plt.show()

