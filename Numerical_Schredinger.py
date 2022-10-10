import numpy as np
import PySimpleGUI as sg
from PIL import ImageGrab


def pit_potential_function(x, pit_height, pit_width, barrier_height=0, barrier_width=0):
    if not (barrier_height and barrier_width):
        return pit_height if abs(x) > pit_width else 0
    else:
        if abs(x) < barrier_width:
            return barrier_height
        elif barrier_width <= abs(x) <= pit_width:
            return 0
        else:
            return pit_height


def hamilton_matrix(interval, N, b, pit_height, pit_width, m, hbar, barrier_height, barrier_width):
    matrix = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            if i == j - 1:
                matrix[i][j] = (-hbar ** 2) / (2 * m * b ** 2) * 7.67
            elif i == j:
                matrix[i][j] = ((hbar ** 2) / (m * b ** 2) * 7.67 + pit_potential_function(i * b + interval[0], pit_height, pit_width, barrier_height, barrier_width))
            elif i == j + 1:
                matrix[i][j] = (-hbar ** 2) / (2 * m * b ** 2) * 7.67
    return matrix


def eigenvalues(hamilton):
    return np.linalg.eigh(hamilton)


def pit_energy_levels(eigenvalues, pit_height):
    index, pit_values = 0, []
    while index != len(eigenvalues) and eigenvalues[index] < pit_height:
        pit_values.append(eigenvalues[index])
        index += 1
    return np.array(pit_values), index


def energy_levels_find(interval, number_of_points, b, pit_height, pit_width, m, hbar, barrier_height, barrier_width):
    matrix = hamilton_matrix(interval, number_of_points, b, pit_height, pit_width, m, hbar, barrier_height, barrier_width)
    
    eigenvalues_matrix, values_matrix = eigenvalues(matrix)
    values_matrix = np.transpose(values_matrix)
    
    pit_energy, stop_index = pit_energy_levels(eigenvalues_matrix, pit_height)
    values_matrix = values_matrix[:stop_index]
    return pit_energy, values_matrix


def draw_x_axis(graph, pit_width, figures_id):
    x_interval = np.linspace(-pit_width * 2, pit_width * 2, 11)
    figures_id['x_axis'] = []
    figures_id['x_axis_txt'] = []
    for index, x in enumerate(range(-200, 201, 40)):
        figures_id['x_axis'].append(graph.DrawLine((x, -203), (x, -197)))
        figures_id['x_axis_txt'].append(graph.DrawText('{0:.1f}'.format(x_interval[index]), (x, -210), color='green'))


def draw_y_axis(graph, pit_height, figures_id):
    y_interval = np.linspace(0, pit_height * 2, 11)
    figures_id['y_axis'] = []
    figures_id['y_axis_txt'] = []
    for index, y in enumerate(range(-200, 201, 40)):
        figures_id['y_axis'].append(graph.DrawLine((-3, y), (3, y)))
        if y != -200:
            figures_id['y_axis_txt'].append(graph.DrawText('{0:.2f}'.format(y_interval[index]), (-15, y), color='blue'))


def draw_border(graph, pit_height=0, pit_width=0, barrier_height=0, barrier_width=0, figures_id=None):
    figures_id['border'] = []
    if pit_height and pit_width and barrier_height and barrier_width:
        x_cord_barrier = barrier_width / pit_width * 100
        y_cord_barrier = barrier_height / pit_height * 200
        figures_id['border'].append(graph.DrawLine((-100, -200), (-100, 0), width=3))
        figures_id['border'].append(graph.DrawLine((100, -200), (100, 0), width=3))
        figures_id['border'].append(graph.DrawLine((-200, 0), (-100, 0), width=3))
        figures_id['border'].append(graph.DrawLine((100, 0), (200, 0), width=3))

        figures_id['border'].append(graph.DrawLine((-100, -200), (-x_cord_barrier, -200), width=3))
        figures_id['border'].append(graph.DrawLine((-x_cord_barrier, -200),
                                                   (-x_cord_barrier, -200 + y_cord_barrier), width=3))
        figures_id['border'].append(graph.DrawLine((-x_cord_barrier, -200 + y_cord_barrier),
                                                   (x_cord_barrier, -200 + y_cord_barrier), width=3))
        figures_id['border'].append(graph.DrawLine((x_cord_barrier, -200 + y_cord_barrier),
                                                   (x_cord_barrier, -200), width=3))
        figures_id['border'].append(graph.DrawLine((x_cord_barrier, -200), (100, -200), width=3))
    else:
        figures_id['border'].append(graph.DrawLine((-100, -200), (-100, 0), width=3))
        figures_id['border'].append(graph.DrawLine((100, -200), (100, 0), width=3))
        figures_id['border'].append(graph.DrawLine((-200, 0), (-100, 0), width=3))
        figures_id['border'].append(graph.DrawLine((100, 0), (200, 0), width=3))
        figures_id['border'].append(graph.DrawLine((-100, -200), (100, -200), width=3))


def delete_graph_lines(element, graph):
    if isinstance(element, list):
        [graph.delete_figure(element[index]) for index in range(len(element))]
    else:
        if element:
            graph.delete_figure(element)


def draw_energy_levels(graph, levels_energy, pit_height, figures_id):
    energy_levels = [-200 + levels_energy[index] / pit_height * 200 for index in range(levels_energy.shape[0])]
    figures_id['energy_levels'] = [graph.DrawLine((-100, energy_levels[index]), (100, energy_levels[index]), 
                                                  color='blue', width=2) for index in range(len(energy_levels))]
    return energy_levels


def wavefunction_array(x_axis, eigenvector):
    return [(x_axis[index], eigenvector[index]) for index in range(len(eigenvector))]


def draw_function(mode, level, graph, figures_id, x_axis, eigenvectors, coord_level, pit_width, scaling, color):
    if mode == 'Wavefunction':
        figures_id['wavefunction'].append(graph.draw_lines(
            wavefunction_array(x_axis * 100 / pit_width,
                               eigenvectors[level] * scaling + coord_level), color=color, width=3))
    else:
        figures_id['wavefunction'].append(graph.draw_lines(
            wavefunction_array(x_axis * 100 / pit_width,
                               np.square(eigenvectors[level]) * scaling + coord_level), color=color, width=3))


def save_element_as_file(element, filename):
    widget = element.Widget

    box = (widget.winfo_rootx(),
           widget.winfo_rooty(),
           widget.winfo_rootx() + widget.winfo_width(),
           widget.winfo_rooty() + widget.winfo_height()
          )

    grab = ImageGrab.grab(bbox=box)
    grab.save(filename)


def redraw_function(figures_id, graph, values, x_axis, eigenvectors, coord_levels, pit_width, scaling_factor, colors):
    delete_graph_lines(figures_id.get('wavefunction', 0), graph)
    if values['-CHOOSEDLEVEL-']:
        draw_levels = list(map(lambda x: int(x) - 1, values['-CHOOSEDLEVEL-']))
        for level in draw_levels:
            draw_function(values['-MODE-'], level, graph, figures_id, x_axis, eigenvectors, coord_levels[level], pit_width, scaling_factor, colors[level])


def check_parameters(window, values, change_configuration_elements):
    check_list = []
    for element in change_configuration_elements:
        if values[element] != '':
            try:
                float(values[element])
                check_list.append(True)
                window[element].update(background_color='White')
            except ValueError:
                check_list.append(False)
                window[element].update(background_color='Yellow')
        else:
            check_list.append(False)
            window[element].update(background_color='Yellow')
    return True if all(check_list) else False


def print_energies(window, pit_energy):
    energies = [f'{pit_energy[index]:.2f} eV\n' for index in range(pit_energy.shape[0])]
    text = ''
    for line in energies:
        text += line
    window['-ENERGYVALUES-'].update('Energy levels values:\n' + text)


def interface():
    number_of_levels = 0
    sg.theme('SystemDefault') # SystemDefaultForReal
    layout = [[sg.Column([[sg.Text('Input pit height: ')],
                          [sg.InputText(default_text=10, border_width=3, size=(15, 1),
                                        key='-PITHEIGHT-'), sg.Text(' eV')],
                          [sg.Text('Input pit width: ')],
                          [sg.InputText(default_text=5, border_width=3, size=(15, 1),
                                        key='-PITWIDTH-'), sg.Text(' A')],
                          [sg.Text('Input mass of particle: ')],
                          [sg.InputText(default_text=1, border_width=3, size=(15, 1),
                                        key='-MASS-'), sg.Text(' * mass of electron')],
                          [sg.Text('For pit with barrier at center: ')],
                          [sg.Text('Input barrier height: ')],
                          [sg.InputText(default_text=0, border_width=3, size=(15, 1),
                                        key='-BARRIERHEIGHT-')],
                          [sg.Text('Input barrier width: ')],
                          [sg.InputText(default_text=0, border_width=3, size=(15, 1),
                                        key='-BARRIERWIDTH-')],
                          [sg.Button('Find levels', key='-FIND-', enable_events=True, border_width=3)],
                          [sg.Text('Number of energy levels = ' + str(number_of_levels), key='-LEVELSNUM-')],
                          [sg.Text('Choose level number: ')],
                          [sg.Combo(['Wavefunction', 'Probability density'], default_value='Wavefunction',
                                    enable_events=True, readonly=True, key='-MODE-')],
                          [sg.Listbox(values=[], select_mode='multiple', enable_events=True,
                                      size=(15, 5), key='-CHOOSEDLEVEL-')]],
                         vertical_alignment='top'),
               sg.Graph(canvas_size=(600, 550), graph_bottom_left=(-215, -265), graph_top_right=(265, 235),
                        background_color='white', key='-GRAPH-'),
               sg.Column([[sg.Text('Function scaling', relief=sg.RELIEF_SOLID, border_width=3, pad=((5, 5), (15, 5)))],
                          [sg.InputText(size=(6, 1), key='-SCALING-'), sg.Button('Ok', enable_events=True,
                                                                                 border_width=3, size=(3,1),
                                                                                 key='-OKSCALE-')],
                          [sg.Text('Colors of functions', relief=sg.RELIEF_SOLID, border_width=3,
                                   pad=((5, 5), (15, 5)))],
                          [sg.Combo(values=[], readonly=True, size=(5, 1), key='-FUNCCOLOR-', enable_events=True),
                           sg.InputText(key='-COLORVALUE-', visible=False, enable_events=True),
                           sg.ColorChooserButton(border_width=3, button_text='         ', key='-COLOR-',
                                                 button_color='Green')],
                          [sg.Checkbox('Show splitting', default=False, key='-SPLITTING-', enable_events=True)],
                          [sg.Button('Save graph', border_width=3, enable_events=True, key='-SAVEGRAPH-')],
                          [sg.Text('Energy levels values: ', key='-ENERGYVALUES-')]], vertical_alignment='top')]]
    
    sg.set_options(dpi_awareness=True)
    window = sg.Window('Numerical solution of the Schrödinger equation', layout, element_justification='left',
                       finalize=True, font=('Times New Roman', 12))
    graph = window['-GRAPH-']
    
    graph.DrawLine((-200, -200), (220, -200), width=2)
    graph.DrawLine((0, -200), (0, 220), width=2)

    graph.DrawLine((-7, 220), (0, 230), width=2)
    graph.DrawLine((0, 230), (7, 220), width=2)
    graph.DrawLine((-7, 220), (7, 220), width=2)
    graph.DrawText('Energy, eV', (35, 225), color='black')

    graph.DrawLine((220, -207), (230, -200), width=2)
    graph.DrawLine((230, -200), (220, -193), width=2)
    graph.DrawLine((220, -207), (220, -193), width=2)
    graph.DrawText('Distance, A', (235, -215), color='black')

    figures_id = dict()
    colors = {}
    
    number_of_points = 1000  # parameter
    pit_height, pit_width, m, hbar, barrier_height, barrier_width = 10, 5, 1, 1, 0, 0  # parameters
    pit_width, barrier_width = pit_width / 2, barrier_width / 2
    interval = [-pit_width - 0.8 * pit_width, pit_width + 0.8 * pit_width]  # parameter
    L = interval[1] - interval[0]
    b = L / (number_of_points + 1)
    x_axis = np.linspace(interval[0], interval[1], number_of_points)
    
    pit_energy, eigenvectors = energy_levels_find(interval, number_of_points, b, pit_height, pit_width, m, hbar, barrier_height, barrier_width)
    coord_levels = []
    
    window['-LEVELSNUM-'].update('Number of energy levels = ' + str(len(pit_energy)))
    
    change_configuration_events = ['-PITHEIGHT-', '-PITWIDTH-', '-MASS-', '-BARRIERHEIGHT-', '-BARRIERWIDTH-']
    
    draw_x_axis(graph, pit_width, figures_id)
    draw_y_axis(graph, pit_height, figures_id)
    draw_border(graph, figures_id=figures_id)
    coord_levels = draw_energy_levels(graph, pit_energy, pit_height, figures_id)
    figures_id['wavefunction'] = []
    figures_id['wavefunction'].append(graph.draw_lines(wavefunction_array(x_axis * 100 / pit_width, eigenvectors[0] * 500 + coord_levels[0]), color='green', width=3))
    print_energies(window, pit_energy)
    colors = {index: 'Green' for index in range(len(pit_energy))}
    window['-SCALING-'].update('500')
    scaling_factor = 500
    window['-CHOOSEDLEVEL-'].update(values=list(range(1, len(coord_levels) + 1)), set_to_index=[0])
    window['-FUNCCOLOR-'].update(values=list(range(1, len(coord_levels) + 1)), set_to_index=[0])
    draw_levels = []
    
    while True:
        event, values = window.read(timeout=20)
        if event == sg.WINDOW_CLOSED:
            break
        if event in '-FIND-':
            calculate_levels = check_parameters(window, values, change_configuration_events)
            if not calculate_levels:
                continue
            pit_height, pit_width, m, barrier_height, barrier_width = \
                list(map(float, [values[element] for element in change_configuration_events]))
            pit_width, barrier_width = pit_width / 2, barrier_width / 2
            if barrier_width < pit_width:
                interval = [-pit_width - 0.8 * pit_width, pit_width + 0.8 * pit_width]
                L = interval[1] - interval[0]
                b = L / (number_of_points + 1)
                x_axis = np.linspace(interval[0], interval[1], number_of_points)
                pit_energy, eigenvectors = energy_levels_find(interval, number_of_points, b, pit_height, pit_width, m, hbar, barrier_height, barrier_width)
                window['-LEVELSNUM-'].update('Number of energy levels = ' + str(len(pit_energy)))
                delete_graph_lines(figures_id.get('x_axis', 0), graph)
                delete_graph_lines(figures_id.get('x_axis_txt', 0), graph)
                delete_graph_lines(figures_id.get('y_axis', 0), graph)
                delete_graph_lines(figures_id.get('y_axis_txt', 0), graph)
                delete_graph_lines(figures_id.get('energy_levels', 0), graph)
                delete_graph_lines(figures_id.get('wavefunction', 0), graph)
                delete_graph_lines(figures_id.get('border', 0), graph)
                draw_x_axis(graph, pit_width, figures_id)
                draw_y_axis(graph, pit_height, figures_id)
                draw_border(graph, pit_height, pit_width, barrier_height, barrier_width, figures_id)
                print_energies(window, np.array([0.0], dtype='float64'))
                window['-CHOOSEDLEVEL-'].update(values=[])
                if pit_energy.shape[0]:
                    coord_levels = draw_energy_levels(graph, pit_energy, pit_height, figures_id)
                    colors = {index: 'Green' for index in range(len(pit_energy))}
                    draw_function(values['-MODE-'], 0, graph, figures_id, x_axis, eigenvectors, coord_levels[0], pit_width, scaling_factor, colors[0])
                    print_energies(window, pit_energy)
                    window['-CHOOSEDLEVEL-'].update(values=list(range(1, len(coord_levels) + 1)), set_to_index=[0])
                    window['-FUNCCOLOR-'].update(values=list(range(1, len(coord_levels) + 1)), set_to_index=[0])
                    window['-COLOR-'].update(button_color=colors[0])
                window['-PITWIDTH-'].update(background_color='White')
                window['-BARRIERWIDTH-'].update(background_color='White')
            else:
                window['-PITWIDTH-'].update(background_color='Yellow')
                window['-BARRIERWIDTH-'].update(background_color='Yellow')
        if event == '-CHOOSEDLEVEL-':
            delete_graph_lines(figures_id.get('wavefunction', 0), graph)
            if values['-CHOOSEDLEVEL-']:
                draw_levels = list(map(lambda x: int(x) - 1, values['-CHOOSEDLEVEL-']))
                delete_graph_lines(figures_id.get('energy_levels', 0), graph)
                delete_graph_lines(figures_id.get('wavefunction', 0), graph)
                if pit_energy.shape[0]:
                    coord_levels = draw_energy_levels(graph, pit_energy, pit_height, figures_id)
                    for level in draw_levels:
                        draw_function(values['-MODE-'], level, graph, figures_id, x_axis, eigenvectors, coord_levels[level], pit_width, scaling_factor, colors[level])
        if event == '-MODE-':
            redraw_function(figures_id, graph, values, x_axis, eigenvectors, coord_levels, pit_width, scaling_factor, colors)
        if event == '-OKSCALE-':
            try:
                scaling_factor = float(values['-SCALING-'])
                if scaling_factor > 0:
                    window['-SCALING-'].update(background_color='White')
                    redraw_function(figures_id, graph, values, x_axis, eigenvectors, coord_levels, pit_width,scaling_factor, colors)
                else:
                    window['-SCALING-'].update(background_color='Yellow')
            except ValueError:
                window['-SCALING-'].update(background_color='Yellow')
        if event == '-FUNCCOLOR-':
            window['-COLOR-'].update(button_color=colors[int(values['-FUNCCOLOR-']) - 1])
        if event == '-COLORVALUE-':
            if values['-COLORVALUE-'] != 'None':
                colors[int(values['-FUNCCOLOR-']) - 1] = values['-COLORVALUE-']
                window['-COLOR-'].update(button_color=colors[int(values['-FUNCCOLOR-']) - 1])
                redraw_function(figures_id, graph, values, x_axis, eigenvectors, coord_levels, pit_width, scaling_factor, colors)
        if event == '-SAVEGRAPH-':
            filename = sg.popup_get_file('Save graph', title='Save', save_as=True, file_types=(('PNG', '*.png'),), no_window=True)
            if filename:
                save_element_as_file(graph, filename)
        if event == '-SPLITTING-':
            if values['-SPLITTING-']:
                if barrier_width and barrier_height:
                    figures_id['splitting'] = []
                    for level in range(0, len(pit_energy), 2):
                        if pit_energy[level] < barrier_height:
                            splitting_values = (np.square(eigenvectors[level + 1]) - np.square(eigenvectors[level])) * scaling_factor
                            coord_level = (coord_levels[level + 1] + coord_levels[level]) / 2
                            figures_id['splitting'].append(graph.draw_lines(
                                    wavefunction_array(x_axis * 100 / pit_width, splitting_values + coord_level),
                                    color='Black', width=3))
            else:
                delete_graph_lines(figures_id.get('splitting', 0), graph)

    window.close()


if __name__ == '__main__':
    interface()
