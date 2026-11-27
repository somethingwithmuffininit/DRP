"""Build the MixBench experiment report from the saved figures and tables below."""

import argparse
from io import BytesIO
from pathlib import Path
from xml.sax.saxutils import escape
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT, TA_RIGHT
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import SimpleDocTemplate, Flowable, Paragraph, Table, TableStyle, Spacer, PageBreak
from pdfrw import PdfReader
from pdfrw.buildxobj import pagexobj
from pdfrw.toreportlab import makerl

SOURCE = Path(__file__).resolve().parent
ROOT = SOURCE.parent
FIGURES = SOURCE / 'figures'
PRESENTATION_FIGURES = ROOT / 'presentation' / 'figures'
DEFAULT_OUTPUT = ROOT / 'output' / 'pdf' / 'mixbench-experiment-report.pdf'
FONTS=Path('/usr/share/fonts/truetype/lato')
for name,file in [('Lato','Lato-Regular'),('LatoBold','Lato-Bold'),('LatoItalic','Lato-Italic')]:
    pdfmetrics.registerFont(TTFont(name,str(FONTS/f'{file}.ttf')))
pdfmetrics.registerFontFamily('Lato',normal='Lato',bold='LatoBold',italic='LatoItalic',boldItalic='LatoBold')
ink, teal = map(colors.HexColor,['#17364E','#087F73'])
muted, rule, pale, wash, rust = map(colors.HexColor,['#486172','#CEDDE1','#EAF3F3','#F5F8FA','#924528'])
W=508
S = {
 'body': ParagraphStyle('body',fontName='Lato',fontSize=11.2,leading=15.1,textColor=ink,spaceAfter=9),
 'title': ParagraphStyle('title',fontName='LatoBold',fontSize=35,leading=39,textColor=ink,spaceAfter=12),
 'subtitle': ParagraphStyle('subtitle',fontName='Lato',fontSize=14,leading=18.5,textColor=ink,spaceAfter=15),
 'h1': ParagraphStyle('h1',fontName='LatoBold',fontSize=28,leading=31,textColor=ink,spaceAfter=11),
 'h2': ParagraphStyle('h2',fontName='LatoBold',fontSize=13,leading=16.2,textColor=ink,spaceBefore=5,spaceAfter=5),
 'kicker': ParagraphStyle('kicker',fontName='LatoBold',fontSize=10,leading=13,textColor=teal,spaceAfter=8),
 'caption': ParagraphStyle('caption',fontName='Lato',fontSize=9.25,leading=12.0,textColor=muted,spaceAfter=10,spaceBefore=4),
 'small': ParagraphStyle('small',fontName='Lato',fontSize=9.7,leading=13.0,textColor=ink,spaceAfter=6),
 'cell': ParagraphStyle('cell',fontName='Lato',fontSize=9.5,leading=12.4,textColor=ink),
 'num': ParagraphStyle('num',fontName='Lato',fontSize=9.5,leading=12.4,textColor=ink,alignment=TA_RIGHT),
 'th': ParagraphStyle('th',fontName='LatoBold',fontSize=9.4,leading=12.0,textColor=colors.white),
 'thnum': ParagraphStyle('thnum',fontName='LatoBold',fontSize=9.4,leading=12.0,textColor=colors.white,alignment=TA_RIGHT),
}
def p(content,style='body'): return Paragraph(content,S[style])
story=[]
def add(t,style='body'): story.append(p(t,style))
def heading(number,title,group):
    if number>1: story.append(PageBreak())
    add(group,'kicker')
    add(title, 'title' if number==1 else 'h1')
def sub(t): add(t,'h2')
def cap(number,t): add(f'<b>Figure {number}.</b> {t}','caption')
def table(rows,widths,nums=(),compact=False):
    data=[[p(str(v),'thnum' if j in nums else 'th') if i==0
           else p(str(v),'num' if j in nums else 'cell') for j,v in enumerate(row)] for i,row in enumerate(rows)]
    out=Table(data,colWidths=widths, hAlign='LEFT', repeatRows=1)
    pad=4 if compact else 5
    out.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),ink),
        ('ROWBACKGROUNDS',(0,1),(-1,-1),(colors.white,wash)),('VALIGN',(0,0),(-1,-1),'TOP'),
        ('LINEBELOW',(0,-1),(-1,-1),.65,rule),('LEFTPADDING',(0,0),(-1,-1),7),
        ('RIGHTPADDING',(0,0),(-1,-1),7),('TOPPADDING',(0,0),(-1,-1),pad),
        ('BOTTOMPADDING',(0,0),(-1,-1),pad)]))
    return out
def tab(rows,widths,nums=(),compact=False):
    story.extend((table(rows,widths,nums,compact),Spacer(1,10)))
def call(t):
    panel=Table([[p(t)]],colWidths=[W],hAlign='LEFT')
    panel.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,-1),pale),('BOX',(0,0),(-1,-1),.5,rule),
        ('LINEBEFORE',(0,0),(0,-1),2,teal),('LEFTPADDING',(0,0),(-1,-1),12),
        ('RIGHTPADDING',(0,0),(-1,-1),12),('TOPPADDING',(0,0),(-1,-1),8),
        ('BOTTOMPADDING',(0,0),(-1,-1),8)]))
    story.extend((panel,Spacer(1,10)))
class Plot(Flowable):
    def __init__(self,path,width=W):
        super().__init__()
        self.xobj=pagexobj(PdfReader(str(path),decompress=False).pages[0])
        x1,y1,x2,y2=map(float,self.xobj.BBox)
        self.origin=x1,y1
        self.scale=width/(x2-x1)
        self.width=width
        self.height=(y2-y1)*self.scale
        self.hAlign='CENTER'
    def draw(self):
        self.canv.saveState()
        self.canv.scale(self.scale,self.scale)
        self.canv.translate(-self.origin[0],-self.origin[1])
        self.canv.doForm(makerl(self.canv,self.xobj))
        self.canv.restoreState()
def fig(name, width=W, *, directory=FIGURES):
    story.append(Plot(directory / f'{name}.pdf', width))
def ref(label,url):
    return f'<link href="{escape(url)}" color="#087F73">{escape(label)}</link>'
def numbers(*items): return '<br/>'.join(items)

heading(1,'How fast do random walks mix?','MIXBENCH EXPERIMENT REPORT')
add('Card shuffles, transition weights, direction memory and shortcuts','subtitle')
add('Start with a deck in a known order. After a few shuffles, some outcomes are still easier to predict than they would be with a perfectly shuffled deck. A walk on a graph raises the same question. If you know where the walker started, how long does that information help you predict where it is now?')
fig('riffle')
cap(1,'Exact total variation (TV) for a 52-card deck in the standard random riffle model. The dashed line marks the main cutoff used in this report. Eight riffles are the first to cross it.')
add('A riffle cuts the deck into two packets and interleaves them without disturbing the order within either packet. Repeating the random move produces a walk through possible deck orders. The calculation counts the probabilities of those orders exactly. At seven riffles, the distance from uniform is 0.3341. After eight it is 0.1672. [1]')
sub('What the network experiments found')
tab([
 ['Experiment','What happened'],
 ['Direction memory on a line','On the 100-position line, every initial heading reached the position cutoff by step <b>141</b>. The ordinary walks took <b>1,898</b> steps.'],
 ['Reweighting friendships','On the karate club graph, optimized weights helped the slowest starts. Metropolis was faster for the median start.'],
 ['Choosing a shortcut','The opposite-side link cut the ring result from <b>195 to 61</b> steps. That shortcut increased the error when estimating visits to half the ring.'],
],[160,348])

heading(2,'How the experiments work','SETUP AND METHODS')
sub('What counts as mixed?')
add('<b>The target.</b> The target gives equal probability to every vertex or deck order. Starting in that distribution leaves it unchanged. TV measures the largest difference between the current and target probabilities of any event, such as being in the left half of a line. At TV 0.25, the discrepancy for every event is at most 25 percentage points.')
add('<b>The measurement.</b> For each selected start, begin with all probability at that state. At each step, each vertex passes on a share of its probability. A vertex holding 40% of the total with a 1/4 chance of taking a link sends 10% of the total along that link. Record the first step when all selected starts are within the cutoff. Most comparisons use 0.25 and every start. The full power grid and the repeated one-link trials use fixed samples of starts. The sampling experiments instead run individual trajectories and measure estimation error.')
sub('Four ways to change the walk')
tab([
 ['Approach','What changes and why it might help'],
 ['Optimized weights','Keep the links and assign the same crossing probability in both directions. Tune those probabilities and the chance of staying put to reduce the slowest lingering imbalance.'],
 ['Direction memory','Store a heading as well as a position. The walker usually moves in that heading, turns at an endpoint and occasionally reverses early. Longer runs take it farther along a line.'],
 ['Clockwise shortcut','Move clockwise with probability 1/2. At the two shortcut endpoints, spend 1/4 on the new jump and 1/4 on staying. Elsewhere, spend 1/2 on staying.'],
 ['One-link selection','Add one missing link. Compare random selection, selection from a vertex with fewest neighbors, the Fiedler rule, and exhaustive search for the largest connectivity gain.'],
],[143,365],compact=True)
sub('Why two different eigenvalues appear')
add('<b>Transition eigenvalues</b> measure how much of a pattern of probability excesses and deficits survives a step. A magnitude of 0.9 leaves 90% of that pattern after one step. Weight optimization keeps the target fixed and shrinks the largest magnitude among the other patterns. [2]')
add('<b>Algebraic connectivity</b> scores the links themselves. Two groups joined by a single bridge are connected, yet the bridge creates a bottleneck. A Laplacian compares values at neighboring vertices. Its second-smallest eigenvalue gives the connectivity score. The Fiedler vector goes with that eigenvalue and gives each vertex a coordinate. The Fiedler rule joins an eligible pair far apart in those coordinates, which can put a new link across a bottleneck. [5, 6]')

heading(3,'Changing the weights on a fixed graph','EXPERIMENTS ON TRANSITION PROBABILITIES')
add('Maximum degree gives each link 1 divided by the largest neighbor count on the graph. Metropolis uses the larger count at that link’s two endpoints. The Fastest Mixing Markov Chain (FMMC) optimizes all weights together. Any remaining probability becomes a chance of staying. [2, 7]')
fig('weights',329)
cap(2,'Five-vertex example. A three-link vertex spends all its probability on moves under maximum degree. Under optimized weights it stays put with probability 1/7. The slowest decay factor drops from 2/3 to 3/7, and both baseline walks take 3 steps to TV 0.25 compared with 2 for FMMC.')
sub('The slowest and typical starts can favor different weights')
add('The karate graph records 78 friendship ties among 34 university karate club members. All three walks use those ties and all 34 starting vertices. On the original graph, optimization helps the slowest start most. Metropolis has the lowest median. [8]')
fig('karate',441)
cap(3,'Steps to TV 0.25 for all 34 starts. Each curve sorts its own starts. The height of the right edge compares the slowest starts, and the middle compares typical starts.')
tab([
 ['On the karate graph','Maximum degree','Metropolis','Optimized'],
 ['Original, slowest start','48','35','<b>26</b>'],
 ['Original, median start','21','<b>15</b>','21.5'],
 ['Fiedler link, slowest start','32','25','<b>17</b>'],
 ['Random link, slowest start','46','34','<b>29</b>'],
],[214,108,85,101],(1,2,3),compact=True)

heading(4,'Direction memory on a line','EXPERIMENTS ON THE STATE OF THE WALK')
add('Picture a walker that flips a coin at every step to choose left or right. Many steps undo the previous move. The <b>lift</b> gives each position two states, one per heading. On a line of <i>n</i> positions a rare move, with probability 1/<i>n</i>, reverses the heading and steps back when possible. The usual move keeps going or turns at the end. To compare positions, add the probabilities for both headings at each position. [3]')
sub('All starts on the same line')
add('We propagated every starting position and, for the lifted walk, both initial headings. This gives both a position crossing and a stricter check that includes heading. The 100-position optimization exceeded the study’s 64-vertex solver cap.','small')
tab([
 ['Positions','Max degree','Metropolis','Optimized','Lifted position','Lifted full'],
 ['10','19','19','19','<b>13</b>','14'],
 ['20','76','76','76','<b>27</b>','30'],
 ['50','474','474','474','<b>70</b>','77'],
 ['100','1,898','1,898','Size cap','<b>141</b>','155'],
],[63,82,81,88,104,90],(0,1,2,3,4,5),compact=True)
add('Steps to TV 0.25. The ordinary walks have one state per position.','caption')
sub('Matched endpoint at seven sizes')
fig('path')
cap(4,'Both walkers start at the same end for each size. The lifted walker starts facing into the line. Both axes are logarithmic. Doubling the length roughly doubles the lifted time and quadruples the maximum-degree time over this range.')
add('At 512 positions, reaching position TV 0.10 took <b>1,227 steps</b> with heading and <b>98,328</b> with maximum degree. The fitted exponents were 1.039 and 2.000. Checks against the reference implementation reproduced the lifted curves on six sizes from 10 to 500 positions.','small')

heading(5,'One shortcut, many outcomes','EXPERIMENTS ON RING PLACEMENT')
add('On a clockwise ring, a concentrated group of walkers can travel around together. A two-way shortcut sends some of that probability to another part of the ring. We fixed one endpoint, tried all 29 eligible partners on a 32-position ring, and included all 32 starting positions. [4]')
fig('ring-pairs')
cap(5,'Left, each dot represents one shortcut. Five placements, including the opposite-side link, took 61 steps versus 195 originally. Right, the shaded area spans six random draws per ring size.')
sub('Link policies on the same clockwise ring')
tab([
 ['Policy','Steps to TV 0.25','Range across draws'],
 ['Original ring','195',''],
 ['Opposite-side link','<b>61</b>',''],
 ['Fiedler link','63',''],
 ['Random link, median of 20','71.5','61 to 151'],
 ['Fewest neighbors, median of 20','64','61 to 172'],
 ['Exhaustive connectivity rule','172',''],
],[272,108,128],(1,2),compact=True)
add('All 29 eligible links tied on algebraic connectivity. The exhaustive rule chose the first link in that tie. The rows include deterministic choices and medians over 20 draws for each randomized rule.','small')
sub('Scaling to 256 positions')
add('At 256 positions, the original took 12,435 steps, the opposite-side link took 3,109, and the six random links ranged from 1,086 to 9,393. Four random links beat the opposite-side link. We repeated the size sweep at four TV cutoffs.')

heading(6,'What happens at a bottleneck?','EXPERIMENTS ON BRIDGES AND TREE LINKS')
add('Extra links can give a walk more ways to cross between groups. They can also change the probability of every move under a weighting rule. We tested both effects by varying the number of bridges between two groups of eight and by adding a single link to a binary tree.')
fig('bottlenecks')
cap(6,'Left, the slowest-start result as the bridge count changes. The two baseline weight rules coincide. Right, steps saved over the unedited tree, averaged over 20 draws or start samples with 16 matched starts per trial.')
add('<b>Two dense groups.</b> With one bridge, optimization took 15 steps versus 28 for either baseline. With eight disjoint bridges, the same comparison was 2 versus 3. In this graph, opening the bottleneck shrank the difference between the weight rules.')
add('<b>Binary tree.</b> Fiedler joined two leaves and improved the crossing from 517 to 501 steps. Random links slowed 15 of 20 trials and the fewest-neighbors rule slowed 12 of 20. Each of the 27 slowdowns raised the highest neighbor count from 3 to 4, so maximum-degree weights fell from 1/3 to 1/4 on every link.')
sub('When do the best-connectivity and fastest-mixing links agree?')
add('We tried every missing link on three graphs and kept all ties for both objectives. The final column shows the original and fastest edited crossing for the slowest start.','small')
tab([
 ['Graph','Links tried','Best connectivity','Fastest mixing','Shared winners','Steps'],
 ['12-position clockwise ring','54','54','6','6','28 to 13'],
 ['Two groups of six','35','25','25','25','17 to 9'],
 ['Karate club','483','1','16','1','48 to 32'],
],[151,58,83,76,67,73],(1,2,3,4,5),compact=True)
add('The two dense groups had the same 25 winners for both objectives. Karate’s best-connectivity edge was among the fastest for the slowest start. Its median-start crossing took three more steps than the best edge under that median metric.','small')

heading(7,'Link selection across graph families','EXPERIMENTS ACROSS GRAPH FAMILIES')
add('We tested a dense random model, a growing model where vertices with more neighbors attract links, and binary trees. In the first sweep, the random models used four graphs per size and the tree used one fixed graph per size. The same graph and maximum-degree transition rule were used before and after a link was added. [9, 10]')
fig('families', directory=PRESENTATION_FIGURES)
cap(7,'All-start results for the small synthetic graph sweep. At 100 vertices the dense model gave a median of 3 steps across all policies. In the growing model, the median fell from 41 to 38 with Fiedler. On the 63-vertex tree, crossings were 120 originally, 148 with the fewest-neighbors draw and 111 with Fiedler.')
sub('Paired repeats with a fixed sample of starts')
add('Ten 100-vertex graphs of each random family and 20 repetitions on one 255-vertex tree gave 180 cases. Each comparison held the same 16 starts within that trial. On the dense graphs all walks crossed in two steps, regardless of the improvement in the connectivity score. The large tree results appear in Figure 6.','small')
tab([
 ['100-vertex model','Link rule','Scaled connectivity gain','Mean steps saved'],
 ['Dense random','Fiedler','0.510','0.0'],
 ['Dense random','Exhaustive','0.537','0.0'],
 ['Preferential attachment','Fiedler','0.825','0.5'],
 ['Preferential attachment','Exhaustive','0.906','0.4'],
],[163,102,121,122],(2,3),compact=True)
sub('Size and setup')
add('A broader sweep measured 405 cases across paths, rings, trees and small pairs of cliques. On paths, the optimized and maximum-degree weights track closely in the decay and mixing measurements. Optimization adds setup time. In the tree and joined-group experiments the value of the same one-link budget varies with size and choice of link. The sweep recorded 12 solver size skips and 20 periodic threshold records.','small')
add('The tree and two-clique size sweeps used eight draws per size for each randomized policy. The repository includes the full size and setup curves.','small')

heading(8,'One shortcut on the Glasgow Subway','EXPERIMENTS WITH REAL NETWORKS')
add('We took the 15 stations from the published Glasgow Subway Outer Circle order. Each station is a vertex and each consecutive pair is a link. A hypothetical connection between Kinning Park and Cowcaddens gives the same edited graph to all four walk rules. These are walks on the station graph with one step per transition. [11]')
gl_fig=Plot(FIGURES / 'glasgow.pdf',225)
gl_tab=table([
 ['Walk','Original','With shortcut'],
 ['Clockwise','43','20'],
 ['Maximum degree','43','<b>16</b>'],
 ['Metropolis','43','21'],
 ['Optimized weights','<b>14</b>','16'],
],[119,60,78],(1,2),compact=True)
pair=Table([[gl_fig,gl_tab]],colWidths=[242,266],hAlign='LEFT')
pair.setStyle(TableStyle([('VALIGN',(0,0),(-1,-1),'MIDDLE'),('LEFTPADDING',(0,0),(-1,-1),0),('RIGHTPADDING',(0,0),(-1,-1),0)]))
story.append(pair)
cap(8,'Glasgow station schematic and the worst-start steps to TV 0.25. The dashed link marks the hypothetical shortcut. We used every station as a start.')
call('<b>Read across a row</b> to see what the new link changed for the same walker. The clockwise walk fell from 43 to 20 steps. The optimized walk rose from 14 to 16.<br/><b>Read down a column</b> to compare weights on the same network. The original ring’s best measured walk used optimized weights. On the edited ring, maximum degree tied it.')
sub('What the other networks add to the comparison')
add('In Glasgow, the topology starts as a simple ring. The karate club has groups of people with uneven numbers of friendship ties. On that graph, maximum degree and optimized weights differ by 22 steps for the slowest start, yet Metropolis leads for the median. The tested random link on the karate graph also moved the optimized walk from 26 to 29 steps. Weight optimization targets long-run decay, while a particular TV cutoff can react differently.')
add('The synthetic rings give a controlled version of the link experiment. On an unedited even ring, maximum degree and Metropolis move every step and alternate between parity classes. Those walks cannot reach TV 0.25. On the 32-position ring, optimized weights can stay put and reach the cutoff in 56 steps. That behavior depends on the move rule and the even ring together.')

heading(9,'What changes as a real network grows?','EXPERIMENTS ON NETWORK SIZE')
sub('Full Western States power grid')
add('The grid has 4,941 vertices and 6,594 physical links. Both walks used the same 32 sampled starts. The sample maximum can be below the all-start maximum. The optimized solve was above the vertex cap. [12]')
tab([
 ['On the full grid','Maximum degree','Metropolis'],
 ['Maximum of 32 sampled starts','28,515','<b>7,721</b>'],
 ['Median of 32 sampled starts','18,797.5','<b>5,154.5</b>'],
],[266,121,121],(1,2))
sub('Small pieces of Glasgow and the grid')
add('We took consecutive 4- to 14-station sections from four Glasgow starting stations. On the grid, six connected samples grew from different roots to 12, 24, 48 and 96 vertices, keeping every original link between selected vertices.')
fig('nested-map', 423, directory=PRESENTATION_FIGURES)
cap(9,'The sequence shows one of the six grid growth orders. Highlighted vertices and links enter at the indicated size. The samples overlap because they come from the same grid.')
fig('nested-clean',447)
cap(10,'Crossings at TV 0.25 over every start in each selected subgraph. The shaded grid ranges show differences across the six realized samples. Glasgow sections of equal size share one shape.')
add('At four Glasgow stations, maximum degree took 2 steps and the lift took 4. A separate nine-station check gave 16 and 12. At 96 grid vertices, maximum degree ranged from <b>289 to 584</b> steps. The Fiedler link lowered the median from 395 to 270. Grid optimization ran through 48 vertices in this sweep.','small')

heading(10,'Faster mixing and estimate accuracy','EXPERIMENTS ON SAMPLING')
add('To estimate the share of visits to the first half of a graph, divide visits to that half by the number of steps. The target is 50%. Successive visits can behave alike, so a long trajectory can contain fewer effectively independent observations than its step count suggests.')
fig('sampling',439)
cap(11,'Root mean squared error for estimating visits to the first half, in percentage points. Left, the first three bars use the karate club and the last two use a 32-position ring. Right, timings for pairs of walks on the same line or ring. The vertical axis on the right is logarithmic.')
sub('Equal numbers of steps')
add('We ran 32 trajectories of 20,000 steps for each walk, starting from the target distribution. On karate, errors were 1.73 percentage points for maximum degree, 2.01 for Metropolis and 1.43 for optimized weights. On the ring the shortcut raised the error from 0.37 to 1.02 percentage points even though it reduced the mixing crossing from 195 to 61. Stationary starts isolate the effect of repeated visits from the effect of where the trajectory begins.')
sub('Equal time to set up and sample')
add('We ran five paired trials for each of 21 configurations. Each trajectory yielded a 10- and 60-second reading, including setup. It also estimated visits to one specified low-degree vertex, with true probability 1/32.','small')
tab([
 ['32-position walk','Half, 10 s','Half, 60 s','One vertex, 10 s','One vertex, 60 s'],
 ['Ordinary line','0.0407','0.0151','0.00495','0.00284'],
 ['Lifted line','<b>0.00688</b>','<b>0.00291</b>','<b>0.000706</b>','<b>0.000435</b>'],
 ['Clockwise ring','<b>0.00447</b>','<b>0.00175</b>','<b>0.00117</b>','<b>0.000324</b>'],
 ['Clockwise ring with shortcut','0.0159','0.00441','0.00167','0.000379'],
],[184,73,73,89,89],(1,2,3,4),compact=True)
add('Error in percentage points. Root mean squared error squares each run’s error, averages those squares, then takes the square root.','caption')
add('The lift helped both line estimates, and the ring shortcut increased both errors. In the original two-group graph at 60 seconds, Metropolis led for left-group visits and optimized weights led for one vertex. These measurements reflect five trials on one machine.','small')

heading(11,'What the comparisons establish','CONCLUSIONS AND SOURCES')
tab([
 ['Question','What the measured comparisons show'],
 ['Can weights alone help?','Yes. On fixed karate ties, optimized weights improved the slowest starts. On a path, the inexpensive baselines already matched the measured optimized crossing.'],
 ['When does direction help?','On longer lines, keeping a heading reduced the tested position-mixing times and the error of both timed estimates. The four-station Glasgow section favored maximum degree.'],
 ['Can connectivity pick the best shortcut?','Sometimes. It identified exactly the fastest links between two cliques. On the ring it tied links with substantially different crossing times.'],
 ['What should a sampling experiment measure?','For a quantity such as the share of visits to half the ring, measure estimation error for that quantity under the intended budget.'],
],[167,341],compact=True)
sub('Experiment inventory')
add('The original audit covered 15 presets and 2,248 planned records, including 105 timed trials. Of those, 2,170 produced measurements and 78 exceeded the solver size limit. Across thresholds and state spaces, 52 records were periodic and zero remained unresolved at the checked horizon. The 180 paired link-policy cases and the 73-case matched ring comparison bring the configured collection to 17 presets. The independent reference calculations supplement those runs.','small')
add('The repository’s <b>EXPERIMENTS.md</b> links the presets and figures. The single <b>notebooks/mixbench.ipynb</b> notebook explains the tests, shows the comparisons and saves the results. It can replay the graph studies, independent reference checks and timed trials. The recorded regression run passed all 16 suites.','small')
sub('Papers and data')
sources=[
 ('1','Bayer and Diaconis (1992). Trailing the Dovetail Shuffle to Its Lair.','https://www.stat.berkeley.edu/~aldous/157/Papers/bayer_diaconis.pdf'),
 ('2','Boyd, Diaconis and Xiao (2004). Fastest Mixing Markov Chain on a Graph.','https://web.stanford.edu/~boyd/papers/fmmc.html'),
 ('3','Diaconis, Holmes and Neal (2000). Analysis of a Nonreversible Markov Chain Sampler.','https://projecteuclid.org/euclid.aoap/1019487508'),
 ('4','Gerencsér (2019). Analysis of a Non-Reversible Markov Chain Speedup by a Single Edge.','https://arxiv.org/abs/1905.03223'),
 ('5','Fiedler (1973). Algebraic Connectivity of Graphs.','https://dml.cz/handle/10338.dmlcz/101168'),
 ('6','Wang and Van Mieghem (2008). Algebraic Connectivity Optimization via Link Addition.','https://www.nas.ewi.tudelft.nl/people/Piet/papers/Bionetics2008_AlgebraicConnectivity.pdf'),
 ('7','Hastings (1970). Monte Carlo Sampling Methods Using Markov Chains and Their Applications.','https://doi.org/10.1093/biomet/57.1.97'),
 ('8','Zachary (1977). An Information Flow Model for Conflict and Fission in Small Groups.','https://doi.org/10.1086/jar.33.4.3629752'),
 ('9','Gilbert (1959). Random Graphs.','https://doi.org/10.1214/aoms/1177706098'),
 ('10','Barabási and Albert (1999). Emergence of Scaling in Random Networks.','https://doi.org/10.1126/science.286.5439.509'),
 ('11','Strathclyde Partnership for Transport. Subway Maps and Stations.','https://www.spt.co.uk/travel-with-spt/subway/maps-stations/'),
 ('12','Watts and Strogatz (1998). Collective Dynamics of Small-World Networks.','https://doi.org/10.1038/30918'),
]
add(numbers(*[f'[{i}] {ref(label,url)}' for i,label,url in sources]),'small')

class Report(SimpleDocTemplate):
    def afterFlowable(self,flowable):
        if isinstance(flowable,Paragraph) and flowable.style.name in ('h1','title'):
            bookmark='report-page-'+str(self.page)
            self.canv.bookmarkPage(bookmark)
            self.canv.addOutlineEntry(flowable.getPlainText(),bookmark,0)

def chrome(canvas,doc):
    w,h=letter
    canvas.setStrokeColor(rule)
    canvas.setLineWidth(.6)
    canvas.line(52,48,w-52,48)
    canvas.setFont('LatoBold',9)
    canvas.setFillColor(teal)
    canvas.drawString(52,31,'MixBench')
    canvas.setFillColor(muted)
    canvas.setFont('Lato',9)
    canvas.drawRightString(w-52,31,str(doc.page))

def build_report(destination):
    doc = Report(
        destination, pagesize=letter, leftMargin=46, rightMargin=46,
        topMargin=45, bottomMargin=64,
        title='How fast do random walks mix? | MixBench experiment report',
        author='MixBench',
        subject='Experiments with card shuffling, Markov chains, graph links and sampling',
    )
    doc.build(story, onFirstPage=chrome, onLaterPages=chrome)


def add_student_credits(source, output):
    from pypdf import PdfWriter
    from pypdf.annotations import FreeText
    from pypdf.generic import (
        ArrayObject, DecodedStreamObject, DictionaryObject, FloatObject, NameObject,
        NumberObject, TextStringObject,
    )

    weights = 'Danielle (weights)'
    direction = 'Dana (direction memory)'
    shortcut = 'Caroline (ring shortcut)'
    links = 'Aditi (link selection)'
    all_students = (weights, direction, shortcut, links)
    credits_by_section = {
        'How fast do random walks mix?': all_students,
        'How the experiments work': all_students,
        'Changing the weights on a fixed graph': (weights,),
        'Direction memory on a line': (direction, weights),
        'One shortcut, many outcomes': (shortcut, links),
        'What happens at a bottleneck?': (links, weights),
        'Link selection across graph families': (links,),
        'One shortcut on the Glasgow Subway': (shortcut, weights),
        'What changes as a real network grows?': (links, direction, weights),
        'Faster mixing and estimate accuracy': (direction, shortcut, weights),
        'What the comparisons establish': all_students,
    }
    writer = PdfWriter(clone_from=source)
    font = writer._add_object(DictionaryObject({
        NameObject('/Type'): NameObject('/Font'),
        NameObject('/Subtype'): NameObject('/Type1'),
        NameObject('/BaseFont'): NameObject('/Helvetica'),
    }))
    for section in writer.outline:
        if section.title not in credits_by_section:
            continue
        page_index = writer.get_destination_page_number(section)
        text = 'Student credit: ' + ' | '.join(credits_by_section[section.title])
        font_size = 8
        width = pdfmetrics.stringWidth(text, 'Helvetica', font_size) + 8
        left = (letter[0] - width) / 2
        note = FreeText(
            text=text, rect=(left, 28, left + width, 40), font_size=f'{font_size}pt',
            font_color='087F73', border_color=None, background_color=None,
        )
        # An explicit appearance keeps the removable annotation visible across PDF viewers.
        appearance = DecodedStreamObject()
        appearance.update({
            NameObject('/Type'): NameObject('/XObject'),
            NameObject('/Subtype'): NameObject('/Form'),
            NameObject('/BBox'): ArrayObject([
                NumberObject(0), NumberObject(0), FloatObject(width), NumberObject(12),
            ]),
            NameObject('/Resources'): DictionaryObject({
                NameObject('/Font'): DictionaryObject({NameObject('/Helv'): font}),
            }),
        })
        escaped = text.replace('\\', '\\\\').replace('(', r'\(').replace(')', r'\)')
        appearance.set_data(
            f'BT /Helv {font_size} Tf 0.0314 0.4980 0.4510 rg 4 3 Td ({escaped}) Tj ET'.encode('ascii')
        )
        note[NameObject('/AP')] = DictionaryObject({NameObject('/N'): writer._add_object(appearance)})
        note[NameObject('/F')] = NumberObject(4)
        note[NameObject('/T')] = TextStringObject('MixBench student credit')
        writer.add_annotation(page_number=page_index, annotation=note)
    writer.write(output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument('--no-student-credits', action='store_true')
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if args.no_student_credits:
        build_report(str(args.output))
    else:
        with BytesIO() as source:
            build_report(source)
            source.seek(0)
            add_student_credits(source, args.output)
