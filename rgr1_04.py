import math
import itertools
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import random # Добавлен импорт для генерации случайных точек
import tkinter as tk # Импорт для графического интерфейса
from tkinter import messagebox # Импорт для стандартных диалоговых окон

# --- Геометрические вспомогательные функции ---
# (Эти функции остаются без изменений, комментарии сохранены из предыдущей версии)

def cross_product(p1, p2, p3):
    """
    Вычисляет векторное произведение (p2-p1) x (p3-p1).
    Знак результата указывает на ориентацию точки p3 относительно вектора p1->p2.
    - Если результат > 0, то p3 находится слева от вектора p1->p2 (поворот против часовой стрелки).
    - Если результат < 0, то p3 находится справа от вектора p1->p2 (поворот по часовой стрелке).
    - Если результат = 0, то точки p1, p2, p3 коллинеарны.
    """
    return (p2[0] - p1[0]) * (p3[1] - p1[1]) - \
           (p2[1] - p1[1]) * (p3[0] - p1[0])

def distance_sq(p1, p2):
    """
    Вычисляет квадрат евклидова расстояния между двумя точками p1 и p2.
    Использование квадрата расстояния позволяет избежать вызова math.sqrt,
    что может быть полезно для сравнения расстояний.
    """
    return (p1[0] - p2[0])**2 + (p1[1] - p2[1])**2

def convex_hull_graham_scan(points):
    """
    Строит минимальную выпуклую оболочку для набора точек с помощью алгоритма сканирования Грэхема.
    """
    if len(points) < 3:
        return sorted(points, key=lambda p: (p[0], p[1]))
    p0 = min(points, key=lambda p: (p[1], p[0]))
    other_points = [p for p in points if p != p0]
    def angle_sort_key(point):
        angle = math.atan2(point[1] - p0[1], point[0] - p0[0])
        return (angle, distance_sq(p0, point))
    sorted_angular_points = sorted(other_points, key=angle_sort_key)
    if not sorted_angular_points:
        return [p0]
    hull = [p0, sorted_angular_points[0]]
    for i in range(1, len(sorted_angular_points)):
        pt = sorted_angular_points[i]
        while len(hull) >= 2 and cross_product(hull[-2], hull[-1], pt) <= 0:
            hull.pop()
        hull.append(pt)
    return hull

def polygon_area(vertices):
    """
    Вычисляет площадь многоугольника по формуле шнурков (Shoelace formula).
    """
    n = len(vertices)
    if n < 3:
        return 0.0
    area = 0.0
    for i in range(n):
        j = (i + 1) % n
        area += vertices[i][0] * vertices[j][1]
        area -= vertices[j][0] * vertices[i][1]
    return abs(area) / 2.0

def is_inside_convex_polygon(point, polygon_vertices):
    """
    Проверяет, находится ли точка `point` внутри или на границе выпуклого многоугольника `polygon_vertices`.
    """
    n = len(polygon_vertices)
    if n < 3:
        return False
    expected_sign = 0 
    epsilon = 1e-9
    for i in range(n):
        p1 = polygon_vertices[i]
        p2 = polygon_vertices[(i + 1) % n]
        cp = cross_product(p1, p2, point)
        current_sign = 0
        if cp > epsilon:
            current_sign = 1
        elif cp < -epsilon:
            current_sign = -1
        if current_sign == 0:
            min_x, max_x = min(p1[0], p2[0]), max(p1[0], p2[0])
            min_y, max_y = min(p1[1], p2[1]), max(p1[1], p2[1])
            if not (min_x - epsilon <= point[0] <= max_x + epsilon and \
                    min_y - epsilon <= point[1] <= max_y + epsilon):
                return False
            continue
        if expected_sign == 0:
            expected_sign = current_sign
        elif expected_sign != current_sign:
            return False
    return True

# --- Функция для генерации случайных точек ---
def generate_random_points(num_points, min_coord, max_coord):
    """Генерирует список случайных точек (x, y)."""
    points = []
    for _ in range(num_points):
        x = random.randint(min_coord, max_coord)
        y = random.randint(min_coord, max_coord)
        points.append((x, y))
    return points

# --- Основная логика поиска четырехугольника ---
def find_min_area_covering_quadrilateral(set1, set2, output_messages_list):
    """
    Находит 4 различные точки из `set1`, образующие выпуклый четырехугольник
    минимальной площади, покрывающий все точки из `set2`.
    """
    min_area = float('inf')
    optimal_quad_vertices = None
    if len(set1) < 4:
        msg = "Предупреждение: В Множестве 1 менее 4 точек. Невозможно сформировать четырехугольник."
        print(msg) # Вывод в консоль
        if msg not in output_messages_list:
             output_messages_list.append(msg)
        return None, float('inf')
    for p_combination in itertools.combinations(set1, 4):
        current_points_list = list(p_combination)
        hull_vertices = convex_hull_graham_scan(current_points_list)
        if len(hull_vertices) != 4:
            continue
        current_area = polygon_area(hull_vertices)
        all_points_in_set2_covered = True
        if not set2:
            pass
        else:
            for point_s2 in set2:
                if not is_inside_convex_polygon(point_s2, hull_vertices):
                    all_points_in_set2_covered = False
                    break
        if all_points_in_set2_covered:
            if current_area < min_area:
                min_area = current_area
                optimal_quad_vertices = hull_vertices
    return optimal_quad_vertices, min_area

# --- Функция для обработки выбора и отображения ---
def process_choice(choice_index, root_window):
    """
    Обрабатывает выбор пользователя, запускает вычисления и отображение.
    """
    global set1_points_global, set2_points_global # Используем глобальные переменные для передачи данных
    
    output_messages = [] # Локальный список сообщений для этого запуска
    
    # Параметры для случайной генерации (можно вынести в константы или сделать настраиваемыми)
    NUM_RANDOM_SET1 = 15
    NUM_RANDOM_SET2 = 5
    COORD_MIN = 0
    COORD_MAX = 100

    # Предопределенные наборы точек
    predefined_sets = [
        { # Пример 1 (основной)
            "name": "Предопределенный 1 (основной)",
            "set1": [(0, 0), (10, 0), (10, 10), (0, 10), (1, 5), (5, 1), (9, 5), (5, 9), (2,2), (8,2), (8,8), (2,8)],
            "set2": [(3, 3), (7, 7), (4, 6), (6, 4)]
        },
        { # Пример 2 (другой)
            "name": "Предопределенный 2 (другой)",
            "set1": [(0,0), (1,3), (4,4), (5,1), (2,-1), (3,0)],
            "set2": [(2,2), (3,1)]
        },
        { # Пример 3 (нет решения)
            "name": "Предопределенный 3 (нет решения)",
            "set1": [(0,0), (1,0), (0,1), (1,1)],
            "set2": [(5,5)]
        },
        { # Пример 4 (мало точек)
            "name": "Предопределенный 4 (мало точек в М1)",
            "set1": [(0,0), (1,0), (0,1)],
            "set2": [(0.5,0.5)]
        }
    ]

    current_set1 = []
    current_set2 = []
    mode_description = ""

    if choice_index == 0: # Случайные точки
        mode_description = f"Режим: Случайные точки ({NUM_RANDOM_SET1} в М1, {NUM_RANDOM_SET2} в М2, коорд. {COORD_MIN}-{COORD_MAX})"
        print(f"\n{mode_description}") # Вывод в консоль
        output_messages.append(mode_description)
        current_set1 = generate_random_points(NUM_RANDOM_SET1, COORD_MIN, COORD_MAX)
        current_set2 = generate_random_points(NUM_RANDOM_SET2, COORD_MIN, COORD_MAX)
        msg1_console = f"Сгенерировано {len(current_set1)} случайных точек для Множества 1."
        msg2_console = f"Сгенерировано {len(current_set2)} случайных точек для Множества 2."
        print(msg1_console) # Вывод в консоль
        print(msg2_console) # Вывод в консоль
        # Эти сообщения не добавляем в messagebox, т.к. mode_description уже информативен
    else: # Предопределенные наборы
        selected_set_data = predefined_sets[choice_index - 1]
        mode_description = f"Режим: {selected_set_data['name']}"
        print(f"\n{mode_description}") # Вывод в консоль
        output_messages.append(mode_description)
        current_set1 = selected_set_data["set1"]
        current_set2 = selected_set_data["set2"]
        print(f"Используется набор: {selected_set_data['name']}") # Вывод в консоль
        print(f"Множество 1: {current_set1}") # Вывод в консоль
        print(f"Множество 2: {current_set2}") # Вывод в консоль


    # Сохраняем выбранные/сгенерированные наборы в глобальные переменные для доступа из Matplotlib
    set1_points_global = current_set1
    set2_points_global = current_set2

    # Вызов основной функции для поиска четырехугольника
    optimal_quad, min_area_val = find_min_area_covering_quadrilateral(current_set1, current_set2, output_messages)

    # Формирование и вывод результатов
    if optimal_quad:
        msg_res1 = f"Найден оптимальный четырехугольник с вершинами: {optimal_quad}"
        msg_res2 = f"Минимальная площадь: {min_area_val:.2f}"
        print(msg_res1) # Вывод в консоль
        print(msg_res2) # Вывод в консоль
        output_messages.append(msg_res1)
        output_messages.append(msg_res2)
    else:
        general_fail_msg = "Не удалось найти четырехугольник, удовлетворяющий условиям."
        # Сообщение о недостатке точек уже добавлено в find_min_area_covering_quadrilateral
        if not any("Предупреждение: В Множестве 1 менее 4 точек" in msg for msg in output_messages) and \
           not any(general_fail_msg in msg for msg in output_messages):
            print(general_fail_msg) # Вывод в консоль
            output_messages.append(general_fail_msg)
    
    final_message_for_box = "\n".join(output_messages)

    # Закрываем окно выбора Tkinter перед показом messagebox и графика
    if root_window:
        root_window.destroy()

    # Отображаем messagebox с результатами
    if final_message_for_box:
        # Создаем временное корневое окно Tkinter для messagebox
        temp_root_for_msgbox = tk.Tk()
        temp_root_for_msgbox.withdraw()  # Скрываем его
        messagebox.showinfo("Результаты обработки", final_message_for_box)
        temp_root_for_msgbox.destroy()

    # --- Графическое представление результатов ---
    if current_set1 or current_set2: # Рисуем, если есть хотя бы одно множество точек
        plt.figure(figsize=(10, 10))
        if current_set1:
            s1_x, s1_y = zip(*current_set1)
            plt.scatter(s1_x, s1_y, c='blue', marker='o', label='Множество 1', s=50, zorder=2)
        if current_set2:
            s2_x, s2_y = zip(*current_set2)
            plt.scatter(s2_x, s2_y, c='red', marker='x', label='Множество 2', s=70, zorder=3)
        if optimal_quad:
            polygon_patch = Polygon(optimal_quad, closed=True, edgecolor='green', facecolor='lightgreen', alpha=0.5, linewidth=2, zorder=1)
            plt.gca().add_patch(polygon_patch)
            for i, (x, y) in enumerate(optimal_quad):
                plt.annotate(f"V{i+1}\n({x},{y})", (x, y), textcoords="offset points", xytext=(5,5), ha='left', fontsize=8)
        plt.xlabel("X координата")
        plt.ylabel("Y координата")
        plt.title("Минимальный покрывающий четырехугольник")
        plt.legend()
        plt.grid(True)
        plt.axis('equal')
        plt.tight_layout()
        plt.show()
    elif final_message_for_box: # Если были сообщения, но нет точек для графика
         print("Нет данных для отображения на графике.")


# --- Создание GUI для выбора ---
def create_choice_window():
    """
    Создает и отображает окно Tkinter для выбора режима.
    """
    root = tk.Tk()
    root.title("Выбор набора данных")
    root.geometry("350x250") # Размер окна

    tk.Label(root, text="Выберите источник данных:").pack(pady=10)

    # Кнопка для случайных точек
    btn_random = tk.Button(root, text="Случайные точки", command=lambda: process_choice(0, root))
    btn_random.pack(pady=5)

    # Кнопки для предопределенных наборов
    predefined_set_names = [
        "Предопределенный 1 (основной)",
        "Предопределенный 2 (другой)",
        "Предопределенный 3 (нет решения)",
        "Предопределенный 4 (мало точек в М1)"
    ]
    for i, name in enumerate(predefined_set_names):
        # Используем lambda с аргументом по умолчанию для передачи правильного индекса
        btn = tk.Button(root, text=name, command=lambda index=i+1: process_choice(index, root))
        btn.pack(pady=2)
    
    root.mainloop()

# --- Точка входа в программу ---
if __name__ == '__main__':
    # Глобальные переменные для хранения выбранных наборов точек,
    # чтобы их можно было использовать в функции отображения Matplotlib,
    # которая вызывается после закрытия окна Tkinter.
    set1_points_global = []
    set2_points_global = []
    
    create_choice_window() # Запускаем GUI для выбора
