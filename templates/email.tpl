JOBID: {{ uuid }} finished.
<p>

{% if has_error %}
There was an error. Please contact us with your job ID for more detail.
{% else %}
Please use the following link to download the resulting pathway. 

<ul>
<li> <a href="{{ url_for('download', uuid=uuid, filename='pathway.pdb', _external=True) }}">pathway.pdb</a> (resulting pathway from the initial and final structures)</li>
<li> <a href="{{ url_for('download', uuid=uuid, filename='close_contacts_5.0_10.0', _external=True) }}">close_contacts_5.0_10.0</a> (residues in close contact in pathway)</li>
</ul>
<p>

To change alignment of pathway please use the following command:

<blockquote>
<pre>
% intra_fit [selection]  (for PyMol)
</pre>
</blockquote>
<p>

Please go to {{ url_for('download', uuid=uuid, _external=True) }} to download the complete output.
{% endif %}

