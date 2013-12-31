JOBID: {{ uuid }} is submitted.
<p>

Please use the following link to check your job status: <a href="{{ url_for('status', uuid=uuid, _external=True) }}">Check Job Status</a>
<p>

When your job is finished, an e-mail will be sent.
<p>
