"""Credit names for GitHub logins (10-17): one place for "who is this".

  credit(st, login) -> 'Name (login)' / 'Name (login, note)' / 'login'

1. site/data/people.yaml (the site's canonical people registry: name, handle, note, aliases): a person whose `handle`
   is the login.  This is where a human says the profile isn't the credit (Kevin Fang for TheSnakeFang, PingYou Ltd
   for BuddenPing, Nate Chaoweeraprasit / SQUISH for itsnaka).
2. the GitHub profile name (users/<login> .name), fetched once and cached in the store's `person` table.
3. the login itself.
The site's tools/credits.py matches people.yaml aliases as substrings, so 'Francisco Couzo (franciscouzo)' resolves.
"""
from __future__ import annotations
import datetime as dt, functools, json, os, subprocess

import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
PEOPLE = os.path.abspath(os.path.join(HERE, '..', '..', '..', 'site', 'data', 'people.yaml'))
SCHEMA = 'CREATE TABLE IF NOT EXISTS person (login TEXT PRIMARY KEY, name TEXT, at TEXT);'


@functools.lru_cache(maxsize=1)
def registry():
    """{login: (name, note)} from people.yaml entries that have a handle"""
    out = {}
    for p in (yaml.safe_load(open(PEOPLE)) or {}).values():
        if isinstance(p, dict) and p.get('handle'):
            out[str(p['handle'])] = (p['name'], p.get('note'))
    return out


def profile_name(st, login):
    """GitHub profile name of login (cached in the store; '' when the profile has none)"""
    st.db.execute(SCHEMA)
    r = st.db.execute('SELECT name FROM person WHERE login=?', (login,)).fetchone()
    if r is not None:
        return r['name']
    from . import github
    try:
        name = json.loads(github._gh([f'users/{login}', '--jq', '{n: .name}']))['n'] or ''
    except (subprocess.CalledProcessError, ValueError, KeyError):
        return ''                                   # not cached: try again next time
    st.db.execute('INSERT OR REPLACE INTO person(login, name, at) VALUES (?,?,?)',
                  (login, name, dt.datetime.now().isoformat(timespec='seconds')))
    return name


def credit(st, login):
    if login in registry():
        name, note = registry()[login]
    else:
        name, note = profile_name(st, login), None
    if not name or name == login:
        return login if not note else f'{login} ({note})'
    return f'{name} ({login}{", " + note if note else ""})'
