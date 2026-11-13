"""Build the Russian PDF report from the reviewed local source changes."""
from pathlib import Path
from xml.sax.saxutils import escape

from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak, Preformatted

OUTPUT = Path(__file__).with_name('psiresp_local_optimization_changes_ru.pdf')
for name, file in [('Body', 'DejaVuSans.ttf'), ('Bold', 'DejaVuSans-Bold.ttf'), ('Mono', 'DejaVuSansMono.ttf')]:
    pdfmetrics.registerFont(TTFont(name, '/usr/share/fonts/truetype/dejavu/' + file))
pdfmetrics.registerFontFamily('Body', normal='Body', bold='Bold', italic='Body', boldItalic='Bold')
ink, blue, muted = colors.HexColor('#172B40'), colors.HexColor('#215B91'), colors.HexColor('#526577')
styles = {
    'title': ParagraphStyle('title', fontName='Bold', fontSize=23, leading=29, textColor=ink, spaceAfter=15),
    'subtitle': ParagraphStyle('subtitle', fontName='Body', fontSize=11, leading=16, textColor=muted, spaceAfter=15),
    'h': ParagraphStyle('h', fontName='Bold', fontSize=14, leading=19, textColor=blue, spaceBefore=13, spaceAfter=8),
    'body': ParagraphStyle('body', fontName='Body', fontSize=10, leading=14, textColor=ink, spaceAfter=7),
    'small': ParagraphStyle('small', fontName='Body', fontSize=8.2, leading=12, textColor=muted, spaceAfter=7),
    'cell': ParagraphStyle('cell', fontName='Body', fontSize=8.6, leading=12, textColor=ink),
    'code': ParagraphStyle('code', fontName='Mono', fontSize=8.1, leading=12, textColor=ink,
                           backColor=colors.HexColor('#F1F5F9'), borderPadding=9, spaceAfter=12),
}
story = []


def p(text, style='body'):
    story.append(Paragraph(text, styles[style]))


def code(text):
    story.append(Preformatted(escape(text), styles['code']))


def table(rows, widths):
    data = [[Paragraph(escape(str(value)), styles['cell']) for value in row] for row in rows]
    item = Table(data, colWidths=widths, repeatRows=1, hAlign='LEFT')
    item.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#DBE8F4')),
        ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#F4F7FA')]),
        ('VALIGN', (0, 0), (-1, -1), 'TOP'),
        ('LEFTPADDING', (0, 0), (-1, -1), 8), ('RIGHTPADDING', (0, 0), (-1, -1), 8),
        ('TOPPADDING', (0, 0), (-1, -1), 7), ('BOTTOMPADDING', (0, 0), (-1, -1), 7),
        ('LINEBELOW', (0, 0), (-1, 0), .6, colors.HexColor('#B8CCDE')),
    ]))
    story.append(item)
    story.append(Spacer(1, 10))


p('Локальная оптимизация<br/>геометрии в psiresp', 'title')
p('Что выполнялось раньше и как работает исправленный путь в ligtools<br/>Состояние локального кода на 13 ноября 2026 года', 'subtitle')
p('Главное изменение', 'h')
p('Одиночный расчёт градиента заменён вызовом <b>psi4.optimize()</b>. Теперь координаты '
  'изменяются итерационно, а результат считается оптимизированным только после успешного '
  'завершения оптимизатора. Изменён локальный путь vendored-форка psiresp 0.4.2 в ligtools.')
table([
    ['Часть расчёта', 'Раньше', 'Теперь'],
    ['Входное задание', 'AtomicInput с driver="gradient"', 'OptimizationInput с исходной молекулой и спецификацией Psi4'],
    ['Выполнение', 'Один расчёт градиента через psi4 --qcschema', 'Скрипт .opt.py вызывает psi4.optimize(); шагами управляет OptKing'],
    ['Координаты результата', 'result.molecule.geometry — исходная геометрия одиночного расчёта', 'result.final_molecule.geometry — геометрия после оптимизации'],
    ['Отметка is_optimized', 'Могла выставляться после одного градиента', 'Выставляется после чтения и проверки результата оптимизации'],
    ['Старые файлы', 'Могли интерпретироваться как оптимизированные структуры', 'Исключены новым ключом кеша'],
], [112, 193, 206])
p('Текущее поведение по умолчанию', 'h')
p('Smi2Gmx использует геометрии MMFF94 без QM-оптимизации. Psi4 рассчитывает волновые функции при фиксированной геометрии. Для включения исправленной QM-оптимизации нужно явно передать optimize_geometry=True при создании Smi2Gmx.')
p('Почему прежний путь был ошибочным', 'h')
p('Градиент показывает направление изменения энергии, но сам по себе не перемещает атомы '
  'и не ищет минимум. В прежнем локальном пути объект AtomicResult возвращал геометрию '
  'одного расчёта; затем существующий Job.optimize_geometries() сохранял её и отмечал '
  'конформер как оптимизированный.')
p('Область исправления', 'h')
p('Исправлен локальный запуск Psi4. Существующий путь QCFractal отдельно отправлял задачу '
  'оптимизации через geometric; его алгоритм этим исправлением не заменялся. '
  'Метод расчёта зарядов остаётся двухстадийным RESP, HF/6-31G* в газовой фазе.')

story.append(PageBreak())
p('Как устроен новый путь', 'title')
p('1. Подготовка задания', 'h')
p('QMGeometryOptimizationOptions.write_input() записывает OptimizationInput: исходную молекулу, '
  'метод, базис и параметры оптимизатора. Создаётся вспомогательный файл .opt.py. Он вызывает '
  'локальный _optimization.py внутри интерпретатора выбранного исполняемого файла Psi4.')
p('2. Настоящая оптимизация', 'h')
code('energy, wavefunction, history = psi4.optimize(\n'
     '    specification["model"]["method"],\n'
     '    molecule=molecule,\n'
     '    return_wfn=True, return_history=True,\n'
     ')')
p('OptKing выбирает перемещения атомов и проверяет сходимость, используя энергии и градиенты '
  'Psi4. Вложенный driver="gradient" сохранён: градиенты нужны на шагах оптимизации. '
  'Он больше не означает, что весь локальный этап ограничен одним градиентом.')
p('3. Сохранение и принятие результата', 'h')
code('final_molecule["geometry"] = np.asarray(\n'
     '    wavefunction.molecule().geometry()\n'
     ')')
p('При успешном возврате psi4.optimize() сохраняются OptimizationResult, конечные координаты, '
  'история энергий и градиентов, число оценок в истории и диагностические величины градиента. '
  'Проверяется отсутствие NaN/Inf. Файл заменяется атомарно через временный файл.')
p('Сходимость проверяет сам оптимизатор; обёртка не реализует независимую повторную проверку '
  'всех критериев. Она ставит маркер converged=True после нормального возврата optimize(). '
  'При исключении записывается success=False и ошибка передаётся вызывающему коду.')
p('При чтении требуется OptimizationResult, маркер локального оптимизатора, подтверждение '
  'сходимости и final_molecule. Одиночный AtomicResult отклоняется. Ошибка процесса Psi4 также '
  'останавливает пайплайн. После принятия результата существующий Job переносит координаты '
  'в конформер и устанавливает is_optimized=True.')
p('Кеш и последующие этапы', 'h')
p('К хешу задания добавлена версия пути psi4-optimize-v1: прежние файлы градиентов не '
  'переиспользуются. Уже созданные в старом ядре Jupyter объекты нужно пересоздать. '
  'run_resp() теперь использует подготовленные ориентации и волновые функции, '
  'не запускает повторно генерацию ориентаций и QM.')

story.append(PageBreak())
p('Файлы, настройки и проверки', 'title')
table([
    ['Файл', 'Что смотреть'],
    ['psiresp/qm.py', 'write_input(), read_output(), _postprocess_result(), get_hash()'],
    ['psiresp/_optimization.py', 'run_optimization(): psi4.optimize(), конечная геометрия, история и ошибки'],
    ['smi2gmx.py', '_run_psi4_input(): запуск .opt.py; run_resp(): использование готовых QM-данных'],
    ['psiresp/job.py', 'optimize_geometries(): существующий перенос результата в конформеры'],
    ['tests/test_local_optimization.py', 'Кеш, отклонение одиночного градиента, сходимость/несходимость и RESP'],
], [189, 322])
p('Настройки библиотеки и настройки прогона TPP', 'h')
table([
    ['Параметр', 'По умолчанию в библиотеке', 'Прогон TPP с QM-оптимизацией'],
    ['Метод / базис', 'HF / 6-31G*', 'HF / 6-31G*'],
    ['max_iter', '200', '200'],
    ['g_convergence (базовый набор)', 'gau_tight', 'gau_tight'],
    ['convergence_tolerance', '1e-3 для четырёх критериев', '1e-3 (последняя попытка)'],
    ['full_hess_every', '10', '−1: модельный гессиан с обновлениями'],
], [143, 184, 184])
p('Значение −1 для TPP было задано только в скрипте прогона после остановки Psi4 системой '
  'из-за нехватки памяти при расчёте полного гессиана. Значение библиотеки осталось 10. '
  'Это изменение настройки оптимизации, а не включение нового протокола расчёта зарядов.')
story.append(PageBreak())
p('Проверки и текущий пайплайн', 'title')
p('На искажённой молекуле воды при строгом gau_tight (convergence_tolerance=None) оптимизация HF/6-31G* прошла 12 шагов: энергия '
  'снизилась с −75,80612412 до −76,01072066 Hartree, конечный максимальный компонент '
  'градиента составил около 7,58 × 10⁻⁸ Hartree/bohr. Проверены изменение координат '
  'и повторное использование сошедшегося кеша. Искусственный лимит в один шаг привёл '
  'к ошибке, без ложного is_optimized=True.')
p('Метанол прошёл всю цепочку в настоящем Jupyter: оптимизация → ESP → RESP → ACPYPE '
  '→ экспорт файлов. Прогон TPP(−2) был остановлен пользователем после 38-го шага '
  '<b>без достижения сходимости</b>; он не подтверждает получение готовых параметров TPP.')
p('При включении QM-оптимизации четыре явных допуска силы и перемещения равны 1e-3 в атомных единицах внутренних координат. convergence_tolerance=None возвращает пороги именованного набора. При линейных координатах это Hartree/bohr и bohr, при угловых — Hartree/рад и рад.')
p('Отдельный полный прогон TPP(−2) на геометрии MMFF94, без QM-оптимизации, успешно прошёл до RTP в копии силового поля: 42 атома, заряд −2. Такой же полный прогон прошли API, MPO и MPI с зарядами −1. Порядок атомов, сохранность геометрии и заряды выходных файлов проверены.')
p('Минимизация RDKit/MMFF94 — отдельный предварительный этап. Она включена по умолчанию '
  'перед окончательным отбором конформеров и не заменяет QM-оптимизацию. '
  'В CoordMatcher оптимизации геометрии нет.', 'small')
p('Основание отчёта: локальные файлы репозитория и журналы ligtools-audit/local-optimization/ '
  'и ligtools-audit/tpp-minus2-full/. Описание относится к рабочей копии, а не к опубликованной версии GitHub.', 'small')


def decorate(canvas, document):
    canvas.saveState()
    canvas.setStrokeColor(colors.HexColor('#CCD9E4'))
    canvas.line(42, 39, A4[0] - 42, 39)
    canvas.setFont('Body', 8)
    canvas.setFillColor(muted)
    canvas.drawString(42, 26, 'ligtools · локальная оптимизация psiresp · 13.11.2026')
    canvas.drawRightString(A4[0] - 42, 26, str(document.page))
    canvas.restoreState()


doc = SimpleDocTemplate(str(OUTPUT), pagesize=A4, leftMargin=42, rightMargin=42,
                        topMargin=39, bottomMargin=53,
                        title='Оптимизация геометрии psiresp: до и после изменений',
                        author='ligtools', subject='Исправление локальной оптимизации Psi4')
doc.build(story, onFirstPage=decorate, onLaterPages=decorate)
print(OUTPUT)
