from ltlf2dfa.Translator import Translator
from ltlf2dfa.DotHandler import DotHandler
import pdb

formula = "G(a->Xb)"
declare_flag = True #True if you want to compute DECLARE assumption for the formula

translator = Translator(formula)
translator.formula_parser()
translator.translate()
translator.createMonafile(declare_flag) #it creates automa.mona file
translator.invoke_mona() #it returns an intermediate automa.dot file

dotHandler = DotHandler()
dotHandler.modify_dot()
dotHandler.output_dot() #it returns the final automa.dot file
# dotHandler.render('automa.dot', view=True)
# 'automa.pdf'
