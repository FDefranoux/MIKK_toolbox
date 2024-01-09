from itertools import product
from collections import deque
import numpy as np
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from bokeh.plotting import figure, save, output_file
from bokeh.models import ColumnDataSource, Plot, Grid, Range1d
from bokeh.models.glyphs import Text, Rect
from bokeh.layouts import gridplot

def needleman_wunsch(x, y):
    """Run the Needleman-Wunsch algorithm on two sequences.

    x, y -- sequences.

    Code based on pseudocode in Section 3 of:

    Naveed, Tahir; Siddiqui, Imitaz Saeed; Ahmed, Shaftab.
    "Parallel Needleman-Wunsch Algorithm for Grid." n.d.
    https://upload.wikimedia.org/wikipedia/en/c/c4/ParallelNeedlemanAlgorithm.pdf
    """
    N, M = len(x), len(y)
    s = lambda a, b: int(a == b)

    DIAG = -1, -1
    LEFT = -1, 0
    UP = 0, -1

    # Create tables F and Ptr
    F = {}
    Ptr = {}

    F[-1, -1] = 0
    for i in range(N):
        F[i, -1] = -i
    for j in range(M):
        F[-1, j] = -j

    option_Ptr = DIAG, LEFT, UP
    for i, j in product(range(N), range(M)):
        option_F = (
            F[i - 1, j - 1] + s(x[i], y[j]),
            F[i - 1, j] - 1,
            F[i, j - 1] - 1,
        )
        F[i, j], Ptr[i, j] = max(zip(option_F, option_Ptr))

    # Work backwards from (N - 1, M - 1) to (0, 0)
    # to find the best alignment.
    alignment = deque()
    i, j = N - 1, M - 1
    while i >= 0 and j >= 0:
        direction = Ptr[i, j]
        if direction == DIAG:
            element = i, j
        elif direction == LEFT:
            element = i, None
        elif direction == UP:
            element = None, j
        alignment.appendleft(element)
        di, dj = direction
        i, j = i + di, j + dj
    while i >= 0:
        alignment.appendleft((i, None))
        i -= 1
    while j >= 0:
        alignment.appendleft((None, j))
        j -= 1

    return list(alignment)

def render_alignment(x, y):
    """Align two sequences, maximizing the
    alignment score, using the Needleman-Wunsch
    algorithm.

    x, y -- sequences.
    """
    alignment = needleman_wunsch(x, y)
    seq_x = "".join(
        "-" if i is None else x[i] for i, _ in alignment
    )
    seq_y = "".join(
        "-" if j is None else y[j] for _, j in alignment
    )
    return '\n'.join([seq_x, seq_y])

def alignements(ref, seq):
    """
    A function to perform alignement between two sequences.
    Args:
        ref (str): reference string
        seq (str): sequence string
        alignement_print (bool): If we want to see the alignement output.

    Returns:
        - score_ratio (float): Ratio between number of hits and len of
        ref sequence
        - prints of formatted alignements
    """
    # Calling funtion for alignements
    alignments = pairwise2.align.globalxx(ref, seq)
    # Calculation of ratio between number of hits and total lenght
    score_ratio = int(alignments[0][2])/max(len(ref), len(seq))
    format = str(format_alignment(*alignments[0]))
    s = format.find('Score')
    print = format[:s]
    return score_ratio, print

def view_alignment(aln, fontsize="9pt", plot_width=800, start=0, ref='HdR_REF'):
    """Bokeh sequence alignment view"""
    #make sequence and id lists from the aln object
    seqs_id = {rec.id: rec.seq for rec in aln}
    seqs = [rec.seq for rec in (aln)]
    ids = [rec.id for rec in aln]    
    text = [i for s in list(seqs) for i in s]
    colors = get_colors(seqs_id, ref=ref)
    N = len(seqs[0])
    S = len(seqs)    
    width = .4

    x = np.arange(1,N+1)
    y = np.arange(0,S,1)
    #creates a 2D grid of coords from the 1D arrays
    xx, yy = np.meshgrid(x, y)
    #flattens the arrays
    gx = xx.ravel()
    gy = yy.flatten()
    #use recty for rect coords with an offset
    recty = gy+.5
    h= 1/S
    #now we can create the ColumnDataSource with all the arrays
    source = ColumnDataSource(dict(x=gx, y=gy, recty=recty, text=text, colors=colors))
    plot_height = len(seqs)*15+50
    x_range = Range1d(0,N+1, bounds='auto')
    if N>100:
        viewlen=100
    else:
        viewlen=N
    #view_range is for the close up view
    view_range = (0,viewlen)
    tools="xpan, xwheel_zoom, reset, save"

    #entire sequence view (no text, with zoom)
    p = figure(title=None, width= plot_width, height=50,
               x_range=x_range, y_range=(0,S), tools=tools,
               min_border=0, toolbar_location='below')
    rects = Rect(x="x", y="recty",  width=1, height=1, fill_color="colors",
                 line_color=None, fill_alpha=0.6)
    p.add_glyph(source, rects)
    p.yaxis.visible = False
    p.grid.visible = False  

    #sequence text view with ability to scroll along x axis
    p1 = figure(title=None, width=plot_width, height=plot_height,
                x_range=view_range, y_range=ids, tools="xpan,reset",
                min_border=0, toolbar_location='below')#, lod_factor=1)          
    glyph = Text(x="x", y="y", text="text", text_align='center',text_color="black",
                text_font_size=fontsize)
    rects = Rect(x="x", y="recty",  width=1, height=1, fill_color="colors",
                line_color=None, fill_alpha=0.4)
    p1.add_glyph(source, glyph)
    p1.add_glyph(source, rects)

    p1.grid.visible = False
    p1.xaxis.major_label_text_font_style = "bold"
    p1.yaxis.minor_tick_line_width = 0
    p1.yaxis.major_tick_line_width = 0
    p1.xaxis.major_label_overrides = {x:str(x+start) for x in range(0,N+1)}
    p.xaxis.major_label_overrides = {x:str(x+start) for x in range(0,N+1)}
    # p1.set_xticklabels(Range1d(start,start+N+1, bounds='auto'), fontsize=20)
    p = gridplot([[p],[p1]], toolbar_location='below')
    return p

def get_colors(seqs, ref='HdR_REF'):
    """Description: determine colors per each variant type from ref in dict of sequences """
    seqs_col = seqs.copy()
    seq_ref = seqs_col.pop(ref)
    dict_colors = {}
    for id, seq in seqs_col.items():
            color_list = []
            assert len(seq) == len(seq_ref), "Won' t work if len are not the same"
            for s in range(0, len(seq)):
                if seq[s] == seq_ref[s]:
                    color_list.append('white')
                elif seq[s] == '-':
                    color_list.append('blue')
                elif seq_ref[s] == '-':
                    color_list.append('red')
                else:
                    color_list.append('green')
            dict_colors[id] = color_list

    ref_colors = {ref:['white' for s in seq_ref]}
    ref_colors.update(dict_colors)
    colors = []
    for l in list(ref_colors.values()):
        colors += l
    return colors


def save_bokeh_figure(figure, filename, title=''):
    """
    Description: Function to save the Bokeh Alignment chart.
    """
    output_file(filename=filename, title=title)
    save(figure)
