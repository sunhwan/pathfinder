JOBID: {{ uuid }} is started.
<p>

Please use the following link to check your job status: <a href="{{ url_for('status', uuid=uuid, _external=True) }}">Check Job Status</a>
<p>

When your job is finished, an e-mail will be sent. Maximum runtime is 12 hours. If you did't get an e-mail after 12 hours, please contact us.

<p>
