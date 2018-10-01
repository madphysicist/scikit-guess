{{fullname | escape | underline }}

.. automodule:: {{fullname}}

{% if functions %}
.. rubric:: Functions

.. autosummary::
    {% for function in functions %}
    {{ function }}
    {% endfor %}

{% for function in functions %}
.. autofunction:: {{ function }}
{% endfor %}

{% endif %}
