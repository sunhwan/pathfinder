JOBID: {{ uuid }} finished.
<p>

{% if has_error %}
There was an error. Please contact us with your job ID for more detail.
{% else %}
Please use the following link to download the resulting pathway. 

<ul>
<li> <a href="{{ url_for('download', uuid=uuid, filename='pathway.pdb', _external=True) }}">pathway.pdb</a> (resulting pathway from the initial and final structures)</li>
<li> <a href="{{ url_for('download', uuid=uuid, filename='movie.gif', _external=True) }}">movie.gif</a> (movie of the resulting pathway)</li>
<li> <a href="{{ url_for('download', uuid=uuid, filename='close_contacts_5.0_10.0', _external=True) }}">close_contacts_5.0_10.0</a> (list of residue paris that were in close contact (< 5 Å) durin the transition but became separated in the end state (> 10 Å))</li>
</ul>
<p>

The transition pathway can be visualized in several molecular viewers, such as PyMOL and VMD. To change alignment of pathway please use the following command:

<blockquote>
<pre>
% intra_fit [selection]  (for PyMol)
</pre>
</blockquote>
<p>

Please go to {{ url_for('download', uuid=uuid, _external=True) }} to download the complete output.
{% endif %}
<p>

Elapsed time:
{% if etime > 3600 %}{{ "%d hours %d mins"|format(etime/3600, (etime-(etime/3600*3600))/60) }}
{% elif etime < 3600 and etime > 60 %}{{ "%d mins"|format(etime/60) }}
{% else %}{{ "%d"|format(etime) }} sec{% endif %}
