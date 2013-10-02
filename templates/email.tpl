JOBID: {{ uuid }} Finished.

{% if has_error %}
There was an error. Please contact us with your job ID for more detail.
{% else %}
Pathway result is attached. To change alignment of pathway please use the following command:

>  % intra_fit [selection]  (for PyMol)

Please go to {{ url_for('download', uuid=uuid, _external=True) }} to download the complete output.
{% endif %}

