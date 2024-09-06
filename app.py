import streamlit as st
from Multiapp import Multipage
import mole,molecule,screen,visual
def apps():
    app = Multipage()
    app.add_page("Molecule Generation and Screening", mole.app)
    app.add_page("Molecule Visualization",visual.app)
    app.add_page("Molecule Selection",molecule.app)
    app.add_page("Molecule Screening",screen.app)
    app.run()
apps()
st.write(" ")
st.write(" ")
st.write(" ")
st.write(" ")
st.write(" ")
st.write(" ")
st.write(" ")
st.write(" ")
st.write(" ")
st.write(" ")
st.write(" ")
st.write(" ")
st.write(" ")



